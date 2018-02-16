pool <- "inst/extdata/pool3pops"

test <- read.pcadapt(pool, type = "pool")
test

tmp <- pcadapt(test)

class(test)
(test + 2)[, 1:5]

test2 <- robust::covRob(t(scale(test, center = TRUE, scale = FALSE))[, -1, drop = FALSE])$dist
plot(test2)
test3 <- pchisq(test2, df = 1.8, lower.tail = FALSE)
test3
plot(test3)
hist(test3)
test4 <- qvalue::qvalue(test3)$qvalue
plot(test4)
hist(test4)


hist(pchisq(test2, df = 1.9, lower.tail = FALSE))

S <- test2
s <- log(seq(exp(sqrt(2)), exp(2), length.out = 100))
df.max <- s[which.max(
  tmp <- sapply(s, function(df) {
    pval <- pchisq(S, df = df, lower.tail = FALSE)
    hist(pval <- pval[pval > 0.2])
    ks.test(pval, "punif", 0.2, 1)$p.value
  }, USE.NAMES = TRUE)
)]

hist(test2, df = df.max, lower.tail = FALSE))

tmp2 <- sapply(s, function(df) {
  median(S) / qchisq(0.5, df = df)
})
plot(tmp2)
