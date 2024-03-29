set.seed(1)
n <- 400
p <- 1000
X <- matrix(rnorm(n*p),n,p)
rownames(X) <- paste0("i",1:n)
colnames(X) <- paste0("j",1:p)
b <- double(p)
b[1:10] <- rnorm(10)
y <- drop(X %*% b + rnorm(n))
fit <- mr_ash(X,y,control = list(max.iter = 500,convtol = 1e-8),
              verbose = "detailed",standardize = TRUE)
elbo <- fit$progress$elbo
elbo <- max(elbo) - elbo + 1e-2
plot(fit$progress$iter,elbo,type = "l",log = "y",lwd = 2,col = "dodgerblue")
plot(b,coef(fit)[-1],pch = 20,col = "darkblue",xlab = "true",ylab = "estimate")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")
