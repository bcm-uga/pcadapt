# ddPCR R package - Dean Attali 2015
# This file contains various UI helper functions for the shiny app

# Create a little question mark link that shows a help popup on hover
helpPopup <- function(content, title = NULL) {
  a(href = "#",
    class = "popover-link",
    `data-toggle` = "popover",
    `data-title` = title,
    `data-content` = content,
    `data-html` = "true",
    `data-trigger` = "hover",
    icon("question-circle")
  )
}
