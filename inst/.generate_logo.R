# +
if (!requireNamespace("hexSticker"))
    install.packages("hexSticker")

if (!requireNamespace("magick"))
    install.packages("magick")
library(magick)

library(hexSticker)
library(plssem)
library(magick)


y <- as.integer(oneIntOrdered$y1)
x <- as.integer(oneIntOrdered$x1)
f <- as.matrix(table(x, y))
d <- f / sum(f)
d


expand_x <- function(x, eps = 1e-5) {
    x0 <- c(x[1] - eps, x[1])
    xm <- x[-c(1)]
    xt <- c(x[length(x)] + 1 - eps, x[length(x)] + 1)

    xm1 <- xm-eps
    xm2 <- xm+eps

    xn <- sort(c(x0, xm1, xm2, xt))
}

expand_table <- function(X, x = seq_len(NROW(X)), y = seq_len(NCOL(X)), eps = 1e-5) {
    xn <- expand_x(x, eps = eps)
    yn <- expand_x(y, eps = eps)

    get_density_single <- function(xi, yi) {
        lowerx <- floor(xi)
        upperx <- ceiling(xi)
        lowery <- floor(yi)
        uppery <- ceiling(yi)

        out <- X[x == lowerx, y == lowery]

        if (!length(out)) 0 else out
    }

    get_density <- function(xn, yn)
        out <- vapply(seq_along(xn), FUN.VALUE=numeric(1L), FUN=\(i) get_density_single(xn[i], yn[i]))

    list(
        z = outer(xn, yn, get_density),
        x = xn,
        y = yn
    )
}

expanded = expand_table(d)

rgbfrom_z <- function(z) {
    d <- z / max(z)# max = 1, min = 0
    rgb(0.85 * d^(1/3) + 0.15, 0.25, 0.5 - 0.7 * z, 0.95)
}

facet_z <- (expanded$z[-1, -1] +
            expanded$z[-nrow(expanded$z), -1] +
            expanded$z[-1, -ncol(expanded$z)] +
            expanded$z[-nrow(expanded$z), -ncol(expanded$z)]) / 4


png(filename = "ord-surface.png", bg = "transparent", width = 6, height = 6, res = 900, units = "in")
persp(x = expanded$x, y = expanded$y, z = expanded$z, theta = 45, phi = 6.2, expand = 0.9, col = rgbfrom_z(facet_z), shade = 0.7, d= 0.9,
     box = F)

dev.off()

sticker("ord-surface.png", package = "", dpi = 900, s_width = 1 + 0.324, s_x = 0.905, s_y = 0.675 + 0.208, h_color = "black", p_size = 35, h_fill = "white",
        p_color = "white", spotlight = FALSE, p_y = 0.5, l_x = 1.3, l_y = 0.7, l_alpha = 0.2, filename = "plssem.png")

plot(image_read("plssem.png"))
