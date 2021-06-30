set.seed(1)
n <- 400
p <- 1000
X <- matrix(rnorm(n*p),n,p)
rownames(X) <- paste0("i",1:n)
colnames(X) <- paste0("j",1:p)
b <- double(p)
b[1:10] <- rnorm(10)
y <- drop(X %*% b + rnorm(n))
fit <- mr.ash(X,y,max.iter = 500,tol = list(convtol = 1e-12))
elbo <- fit$varobj
elbo <- elbo - min(elbo) + 1e-6
plot(1:fit$iter,elbo,type = "l",log = "y",lwd = 2,col = "dodgerblue")
plot(b,coef(fit),pch = 20,col = "darkblue",xlab = "true",ylab = "estimate")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")
