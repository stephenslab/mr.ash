context("init")

test_that(
  "init_mr_ash and glmnet give the same output with the same parameters",
  {
    
    set.seed(7)
    dat <- simulate_regression_data(n = 400, p = 100, s = 4, center = FALSE)
    X <- dat$X
    y <- dat$y
    # Try all four combinations of intercept and standardize.
    for (intercept in c(TRUE, FALSE)) {
      
      for (standardize in c(TRUE, FALSE)) {
        
        set.seed(1)
        fit0 <- init_mr_ash(X, y, init.method = "glmnet", intercept = intercept,
                            standardize = standardize)
        
        mr_ash_coefs <- with(fit0, c(b0, b))
        
        set.seed(1)
        fit.glmnet <- glmnet::cv.glmnet(X, y, intercept = intercept, standardize = standardize)
        
        glmnet_coefs <- as.numeric(coef(fit.glmnet, s = "lambda.1se"))
        
        expect_equal(glmnet_coefs, unname(mr_ash_coefs), tolerance = .15,
                     label = glue::glue("glment coefficients with standardize = {standardize}",
                                        " and intercept = {intercept}"
                                        ),
                     expected.label = glue::glue("mr_ash coefficients with standardize = {standardize}",
                                                 " and intercept = {intercept}"
                     )
                    )
        
      }
      
    }
    
  }
)