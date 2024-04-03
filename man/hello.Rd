\name{easy_table}
\alias{Easy Table}
\title{Easy Table}
\usage{
easy_table()
}
\description{
Prints a regression table as a flextable object.
}
\examples{
test <- lm(mpg ~ factor(cyl), data = mtcars)
test2 <- lm(mpg ~ factor(cyl) + log(hp) + wt + am, data = mtcars)
test3 <- lm(mpg ~ factor(cyl) + log(hp) + wt + log(disp) + am, data = mtcars)
test4 <- glm(am ~ factor(cyl) + log(hp), data = mtcars, family = "binomial")

mymods <- list(
  a1 = test,
  a2 = test2,
  a3 = test3,
  a4 = test4
)

easy_table(mymods,
      highlight = F,
      csv = "test",
      control.var = c("hp", "cyl", "am"),
      robust.se = T
      ) |>
  flextable::width(j = 1, 5, "cm")
}
