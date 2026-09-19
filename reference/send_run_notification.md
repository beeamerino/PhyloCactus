# Send a Run Notification by Email

Sends a short message when a long analysis ends, so that a run measured
in hours does not have to be watched. Intended for the steps this
package runs on a local machine, not on a cluster, where the scheduler
already mails its own notifications.

## Usage

``` r
send_run_notification(
  subject,
  body,
  to = NULL,
  credentials = NULL,
  enabled = TRUE
)
```

## Arguments

- subject:

  Character. Subject line.

- body:

  Character. Message body, interpreted as Markdown by
  [`blastula::md()`](https://rstudio.github.io/blastula/reference/md.html).

- to:

  Character or `NULL`. Recipient, also used as the sender. Defaults to
  `NULL`, which reads the `MY_EMAIL` environment variable, normally set
  in `~/.Renviron`.

- credentials:

  Character or `NULL`. Path to a `blastula` credentials file. Defaults
  to `NULL`, which reads the `PHYLOCACTUS_SMTP_CREDS` environment
  variable and falls back to `~/.blastula_gmail`.

- enabled:

  Logical. `FALSE` returns immediately and sends nothing, so that a
  script can carry the call unconditionally and switch it from one
  place. Defaults to `TRUE`.

## Value

Invisibly `TRUE` when the message was handed to the SMTP server, `FALSE`
otherwise.

## Details

**This function never fails a run.** Every reason it might not send, a
missing address, a missing package, a missing credentials file, an SMTP
error, is reported with a
[`message()`](https://rdrr.io/r/base/message.html) and returns `FALSE`.
An analysis that took sixteen hours must not be lost because a mail
server refused a connection at the end of it.

**Credentials are never read by this package.** The path is handed to
[`blastula::creds_file()`](https://rstudio.github.io/blastula/reference/credential_helpers.html),
which opens it. Nothing from that file is read, printed, logged or
written anywhere by `PhyloCactus`, and the file belongs outside the
repository. Create it once with
[`blastula::create_smtp_creds_file()`](https://rstudio.github.io/blastula/reference/create_smtp_creds_file.html),
and for a Gmail account use an application password rather than the
account password.

## Examples

``` r
if (FALSE) { # \dontrun{
send_run_notification(
  subject = "PhyloCactus finished: treePL",
  body = "The dating run finished. Outputs are in 8_Dating."
)
} # }
```
