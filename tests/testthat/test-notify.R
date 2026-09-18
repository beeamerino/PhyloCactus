# The notification is a convenience attached to a run measured in hours. Two properties matter more
# than what the message says: it must never fail the run, and it must never touch the credentials
# file. Both are tested here. What the message says is tested on the composer, which needs no
# network, no account and no credentials, and is therefore the part that would otherwise never be
# exercised at all.

test_that(".format_elapsed reads as a duration at every scale", {
  expect_equal(.format_elapsed(0), "0s")
  expect_equal(.format_elapsed(45), "45s")
  expect_equal(.format_elapsed(90), "1m 30s")
  expect_equal(.format_elapsed(3600), "1h 00m 00s")
  expect_equal(.format_elapsed(58230), "16h 10m 30s")
  # A negative difference is a clock that moved, not a run that went backwards.
  expect_equal(.format_elapsed(-5), "0s")
})

test_that(".compose_run_notification() reports a finished run with its times and its outputs", {
  started <- as.POSIXct("2026-09-18 01:00:00", tz = "UTC")
  finished <- as.POSIXct("2026-09-18 17:10:30", tz = "UTC")

  note <- .compose_run_notification(
    analysis = "treePL divergence time estimation",
    status = "finished", started = started, finished = finished,
    outputs = c("Best maximum-likelihood chronogram" = "8_Dating/BestTree_treePL.tree",
                "Dating output directory" = "8_Dating"),
    host = "test-host"
  )

  expect_match(note$subject, "^PhyloCactus finished: treePL")
  expect_match(note$body, "16h 10m 30s", fixed = TRUE)
  expect_match(note$body, "8_Dating/BestTree_treePL.tree", fixed = TRUE)
  expect_match(note$body, "Best maximum-likelihood chronogram", fixed = TRUE)
  expect_match(note$body, "test-host", fixed = TRUE)
  # Nothing to report as an error on a run that finished.
  expect_false(grepl("**Error**", note$body, fixed = TRUE))
})

test_that(".compose_run_notification() carries the R error when the run failed", {
  started <- as.POSIXct("2026-09-18 01:00:00", tz = "UTC")

  note <- .compose_run_notification(
    analysis = "treePL divergence time estimation",
    status = "failed", started = started, finished = started + 12,
    outputs = c("Dating output directory" = "8_Dating"),
    error_message = "treePL exited with status 1 for label BS_17",
    host = "test-host"
  )

  # The subject has to be readable on a phone's lock screen, so the failure is in it.
  expect_match(note$subject, "FAILED", fixed = TRUE)
  expect_match(note$body, "treePL exited with status 1 for label BS_17", fixed = TRUE)
  expect_match(note$body, "12s", fixed = TRUE)
})

test_that("send_run_notification() refuses quietly rather than failing the run", {
  # enabled = FALSE is the switch a script leaves in place and turns off from one line.
  expect_silent(res <- send_run_notification("s", "b", to = "someone@example.org",
                                             enabled = FALSE))
  expect_false(res)

  # No recipient is a configuration gap, not an error.
  withr::local_envvar(MY_EMAIL = "")
  expect_message(res <- send_run_notification("s", "b"), "No recipient address")
  expect_false(res)
})

test_that("send_run_notification() reports a missing credentials file and does not stop", {
  skip_if_not_installed("blastula")

  tmp <- withr::local_tempdir()
  missing <- file.path(tmp, "no_such_creds")

  expect_message(
    res <- send_run_notification("s", "b", to = "someone@example.org", credentials = missing),
    "No blastula credentials file"
  )
  expect_false(res)
  # The path is named so the user can act on it, and nothing was created in its place.
  expect_false(file.exists(missing))
})

test_that("automate_treePL() still fails with its own error when notify is on", {
  # The wrapper catches the error to report it and must re-raise it unchanged. A notification that
  # swallowed a failure would be worse than no notification: the run would look successful.
  tmp <- withr::local_tempdir()
  withr::local_envvar(MY_EMAIL = "")

  expect_error(
    suppressMessages(automate_treePL(
      cfg_file = file.path(tmp, "absent.cfg"),
      ml_tree_file = file.path(tmp, "absent.tree"),
      bs_trees_file = file.path(tmp, "absent_bs.tree"),
      results_dir = file.path(tmp, "res"),
      treePL_out = file.path(tmp, "out"),
      notify = TRUE
    ))
  )
})
