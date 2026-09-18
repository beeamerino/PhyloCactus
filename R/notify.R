#' Format an Elapsed Time in Hours, Minutes and Seconds
#'
#' @param secs Numeric. Elapsed seconds.
#' @return A single character string.
#' @noRd
.format_elapsed <- function(secs) {
  secs <- max(0, round(as.numeric(secs)))
  h <- secs %/% 3600
  m <- (secs %% 3600) %/% 60
  s <- secs %% 60
  if (h > 0) {
    sprintf("%dh %02dm %02ds", h, m, s)
  } else if (m > 0) {
    sprintf("%dm %02ds", m, s)
  } else {
    sprintf("%ds", s)
  }
}

#' Compose the Subject and Body of a Run Notification
#'
#' Internal helper. Kept separate from the sending so that what the message says can be tested
#' without a network, an SMTP account or a credentials file, which is the part that would otherwise
#' never be exercised.
#'
#' @param analysis Character. What ran, in words, for the subject line.
#' @param status `"finished"` or `"failed"`.
#' @param started,finished Objects of class `POSIXct` bounding the run.
#' @param outputs Named character vector of paths worth naming in the message. Names become labels.
#' @param error_message Character or `NULL`. The R error, when the run failed.
#' @param host Character. Machine name.
#' @return A list with `subject` and `body`, the latter in Markdown.
#' @noRd
.compose_run_notification <- function(analysis,
                                      status = c("finished", "failed"),
                                      started,
                                      finished = Sys.time(),
                                      outputs = character(0),
                                      error_message = NULL,
                                      host = Sys.info()[["nodename"]]) {
  status <- match.arg(status)

  version <- tryCatch(as.character(utils::packageVersion("PhyloCactus")),
                      error = function(e) "unknown")

  subject <- paste0("PhyloCactus ",
                    if (status == "finished") "finished" else "FAILED",
                    ": ", analysis)

  lines <- c(
    paste0("**", analysis, "** ", status, " on `", host, "`."),
    "",
    paste0("- Started: ", format(started, "%Y-%m-%d %H:%M:%S")),
    paste0("- Finished: ", format(finished, "%Y-%m-%d %H:%M:%S")),
    paste0("- Elapsed: ", .format_elapsed(difftime(finished, started, units = "secs"))),
    paste0("- Host: ", host),
    paste0("- R: ", R.version.string),
    paste0("- PhyloCactus: ", version)
  )

  if (length(outputs) > 0L) {
    labels <- names(outputs)
    if (is.null(labels)) labels <- rep("", length(outputs))
    labels <- ifelse(nzchar(labels), paste0(labels, ": "), "")
    lines <- c(lines, "", "**Outputs**", "",
               paste0("- ", labels, "`", unname(outputs), "`"))
  }

  if (!is.null(error_message) && nzchar(error_message)) {
    lines <- c(lines, "", "**Error**", "", paste0("    ", error_message))
  }

  list(subject = subject, body = paste(lines, collapse = "\n"))
}

#' Send a Run Notification by Email
#'
#' Sends a short message when a long analysis ends, so that a run measured in hours does not have to
#' be watched. Intended for the steps this package runs on a local machine rather than on a cluster,
#' where the scheduler already mails its own notifications.
#'
#' **This function never fails a run.** Every reason it might not send, a missing address, a missing
#' package, a missing credentials file, an SMTP error, is reported with a `message()` and returns
#' `FALSE`. An analysis that took sixteen hours must not be lost because a mail server refused a
#' connection at the end of it.
#'
#' **Credentials are never read by this package.** The path is handed to `blastula::creds_file()`,
#' which opens it. Nothing from that file is read, printed, logged or written anywhere by
#' `PhyloCactus`, and the file belongs outside the repository. Create it once with
#' `blastula::create_smtp_creds_file()`, and for a Gmail account use an application password rather
#' than the account password.
#'
#' @param subject Character. Subject line.
#' @param body Character. Message body, interpreted as Markdown by `blastula::md()`.
#' @param to Character or `NULL`. Recipient, also used as the sender. Defaults to `NULL`, which
#'   reads the `MY_EMAIL` environment variable, normally set in `~/.Renviron`.
#' @param credentials Character or `NULL`. Path to a `blastula` credentials file. Defaults to
#'   `NULL`, which reads the `PHYLOCACTUS_SMTP_CREDS` environment variable and falls back to
#'   `~/.blastula_gmail`.
#' @param enabled Logical. `FALSE` returns immediately and sends nothing, so that a script can carry
#'   the call unconditionally and switch it from one place. Defaults to `TRUE`.
#' @return Invisibly `TRUE` when the message was handed to the SMTP server, `FALSE` otherwise.
#' @examples
#' \dontrun{
#' send_run_notification(
#'   subject = "PhyloCactus finished: treePL",
#'   body = "The dating run finished. Outputs are in 8_Dating."
#' )
#' }
#' @export
send_run_notification <- function(subject, body, to = NULL, credentials = NULL, enabled = TRUE) {
  if (!isTRUE(enabled)) return(invisible(FALSE))

  if (is.null(to)) to <- Sys.getenv("MY_EMAIL", "")
  if (!is.character(to) || length(to) != 1L || !nzchar(to)) {
    message("No recipient address: set MY_EMAIL in ~/.Renviron or pass `to`. No email sent.")
    return(invisible(FALSE))
  }

  if (is.null(credentials)) {
    credentials <- Sys.getenv("PHYLOCACTUS_SMTP_CREDS", "~/.blastula_gmail")
  }
  credentials <- path.expand(credentials)

  if (!requireNamespace("blastula", quietly = TRUE)) {
    message("Package 'blastula' is not installed. No email sent.")
    return(invisible(FALSE))
  }

  if (!file.exists(credentials)) {
    message("No blastula credentials file at '", credentials,
            "'. Create one with blastula::create_smtp_creds_file(). No email sent.")
    return(invisible(FALSE))
  }

  sent <- tryCatch({
    blastula::smtp_send(
      email = blastula::compose_email(body = blastula::md(body)),
      from = to,
      to = to,
      subject = subject,
      credentials = blastula::creds_file(credentials)
    )
    TRUE
  }, error = function(e) {
    message("Email could not be sent, and the run is unaffected: ", conditionMessage(e))
    FALSE
  })

  if (isTRUE(sent)) message("Notification sent: ", subject)
  invisible(isTRUE(sent))
}
