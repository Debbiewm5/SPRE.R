#' @title SPRE: A Semiparametric Data Estimation model to predict outcomes for single-subject and small group designs
#'
#' @description SPRE is a Semiparanetric Data Estimation model to predict outcomes for single-subject and small group designs where the data must be named: Dataset <- read.csv("~/Documents/NameOfYourFile.csv"). There are 2 Excel data columns: FData  Session. SPRE predicts residuals, normal Q-Q, predictions from a change point, effective clinical stability, SPRE Weibull probability, computations for the highest F statistic,
#' ram Session A numberic vector
#' @return  residuals, normal Q-Q, predictions from a change point, effective clinical stability, SPRE Weibull probability. Computations for the highest F statistic, all predictions, the probability curve for effective clinical stability.
#'
#' @example Dataset <- read.csv("~/Documents/HEAT_LILLY_Heat.Stat.14_2016.csv") #' FData 4.8462,4.6429,4.625,4.4375,4.3792,4.375,4.1765,4.1429,4.1333,4.0667,4,4,3.9375,3.875   #' Session 4,3,14,11,13,6,10,8,12,2,1,5,7,9
#'library
#' @depends:
library(MASS)
               save(list=ls(all=TRUE), file="all.Rdata")

#' @name linearDist
#' @description provides a least squares regression (OLS) for the entire dataset from which the coefficient for the slope of the OLS line is obtained for the entire model
#' @param FData A numeric vector
#' @param Session A numeric vector
#' @return beta coefficient A number

    linearDist <- function(lm) {
    lm(Dataset$FData~Dataset$Session)
}
        Session <- Dataset$Session
            stat <- Dataset$FData
            SE <- Dataset$Session
            D <- Dataset$FData
                x <- summary(lm(FData ~ Session, data=Dataset))
                        sL <- x$coefficients[2,1]
                        sL <- abs(sL)

#' @name fDist
#' @description The generic matrix function performs a linear model backwards stepwise regression reducing the model by one row each iteration to determine the highest F statistic in the entire dataset. This value is the 'change point' which indicates where the patient(s)/client(s) adapts to treatment.
#' @param FData A numeric vector
#' @param Session A numeric vector
#' @return fDist[i,]
#' @export

    fDist <- matrix (
        data=NA,
        nrow= length(D),
        ncol=6,
        byrow=TRUE)
            dimnames(fDist) = list(c(rownames.force = NULL),
                    (c("Rsquare","Fstatistic","FstatNumer.",
                    "FstatDenom.","Pvalue","Session")))

            for(i in length(D):1) {
                    x <- summary(lm(FData ~ Session, data=Dataset))
                        Rsq <- x$r.squared
                        n <- x$fstatistic
                        f <- as.numeric(x$fstatistic)
                        P <- as.numeric(pf(f[1],f[2],f[3],lower.tail=F))
                            SE <- Dataset$Session
                            SE <- t(SE)
                                Dataset <- Dataset[-c(nrow(Dataset)),]
                        fDist[i,]= (c(Rsq,n,P,SE[i]))
        }

cat("     Backwards stepwise OLS linear regression ")
cat("\n","\n")
print(signif(fDist,digits=6))
cat(" \n ")

#' @name tester
#' @description Point estimation data from the 'fDist' matrix evaluation
#' @param qn A number
#' @param pv A number
#' @return computations

        tester <- apply(fDist, FUN=which.max,2)
                y <- tester[2]
                    q <- fDist[y,2]
                    qn <- fDist[y,6]
                    v <- fDist[y,6]
                    u <- fDist[y,6]+1
                    w <- fDist[y,1]
                    pv <- fDist[y,5]
                        k <- abs(log(abs(1-(sL*v))))
                saveRDS(qn, "qn.rds")
                saveRDS(pv, "pv.rds")

cat(" \n ")
cat("slope parameter: ")
cat("betaHat1 = "); print(signif(sL,digits=4))
cat(" \n ")
cat("shape parameter: ")
cat("k_value -"); print(signif(k,digits=4))
cat(" \n ")
cat("scale parameter is the value (tau) at the change point - given as the highest F statistic")
cat("\n"); cat("\n")
cat("Session number for - "); print(signif(qn,digits=4))
cat(" \n ")
cat("Session value for - "); print(signif(pv,digits=4))
cat(" \n ")
cat("If R^2 is less than 0.45 at the change point, do not predict for inference")
cat(" \n ")

#' Point estimations from the full saved original dataset

        load("all.Rdata")
            Dataset <- as.data.frame(Dataset)
                z <- Dataset$FData[y]
                sess <- v
        r <- ((1-exp(-(u/v)^k))/(1-exp(-(v/v)^k))) == ((1-exp(-(u/v)^k))/(1-exp(-(v/v)^k)))

#' @name Pred
#' @description The generic matrix function creates the point estimations, which are predictions in this model. The point estimations start at the session after the 'change point'
#' @param fDist A function
#' @param z[i]  A numeric vector
#' @return z[i,]
#' @export

        Pred <- matrix(
                data=NA,
                nrow=40,
                ncol=5,
                byrow=TRUE)
                dimnames(Pred) = list(c(rownames.force = NULL),
                      (c("Session","Ratio","tau","k.factor","Predictions")))

                for (i in 2:40) {
                        if(stat[2] < stat[length(stat)]){
                                v <- i+fDist[y,6]
                                u <- i+(fDist[y,6]+1)
                                sess <- i+fDist[y,6]
                                z[i] <- Dataset$FData[y]*r
                } else {
                        if(stat[2] > stat[length(stat)])
                                v <- i+(fDist[y,6]+1)
                                u <- i+fDist[y,6]
                                sess <- i+fDist[y,6]
                                z[i] <- Dataset$FData[y]*r
        }
                r <-((1-exp(-(u/v)^k))/(1-exp(-(v/v)^k)))
                        z[i] <- z[i-1]*r

                Pred[i,] = c(sess,r,fDist[y,6],k,z[i])
}
cat("\n", "\n")
print(signif(Pred,digits=6))
cat("\n")
cat("Effective clinical stability may be determined when the second significant digit (after the decimal place) remains constant for 2 predictions.
These 2 sessions, or more, are then effective clinical stability as long as treatment is maintained to these sessions.")
cat("\n", "\n")

#' @name myPlots
#' @description The function creates the plots for SPRE residuals, Normal Q-Q and SPRE predictions from a 'change point'
#' @param Y A numeric vector
#' @param S A numeric vector
#' @return plots
#' @export

                myPlots <- function (data,d) {
                        data <-list('Pred',header=TRUE)
                        write.table(data,"data.dat")
                        d <- read.table('data.dat',header=T)
}
                for(i in 1:40){
                        data <- data.frame(data=Pred,
                           nrow=40,
                           ncol=5,
                           byrow=TRUE)
                dimnames(Pred) = list(c(rownames.force = NULL),
                              (c("Session","Ratio","tau","k.factor","Predictions")))
}
                Y <- data$data.Predictions
                S <- data$data.Session
                    lmFit <- lm(Dataset$FData~Dataset$Session, data = Dataset)
        
        dev.new()
        par(mfrow = c(3,1), mar=c(4.7,5.8,2.5,3.8))
        
                plot(lmFit,1,ann=FALSE, pch=18,
                title("SPRE Residuals",col.main="blue",xlab="Fitted Values",
                ylab="Residuals"))
                        abline(0,0)

                plot(lmFit,2,ann=FALSE, pch=18,
                        title("Normal Q-Q",col.main="blue",xlab="xlab
                              Quantiles",
                        ylab="Standardized Residuals"))
                                abline(0,0)

                plot(Y~S, ann=FALSE, type="o", pch=19, col="blue")
                        xlim <- S[-(nrow=1)]
                        ylim <- Y[-(nrow=1)]
                        lines(Y~S,lwd=2, col="blue")
                legend("bottomleft", legend = "ChangePt       ", cex = 0.8)
                title("SPRE Predictions from a 'change point'",col.main="blue",
                        xlab="Prediction Sessions",ylab="Predictions")
                dev.set(dev.prev())
                dev.set(dev.hold())
                
#'  The generic matrix function 'Pred' computes the the second significant digit (after the decimal place) which remains constant for 2 predictions. This is clinical stability.

#' @name Pred
#' @description The generic matrix function creates the point estimations, which are predictions in this model. The point estimations start at the session after the 'change point' to determine the first duplicate predicted values. The associated session is then the effective clinical stability.
#' @param fDist A function
#' @param dup A number
#' @return xu2
#' @export

        Pred <- matrix(
                data=NA,
                nrow=40,
                ncol=1,
                 byrow=TRUE)
                dimnames(Pred) = list(c(rownames.force = NULL),
                      (c("Predictions")))

                for (i in 2:40) {
                    if(stat[2] < stat[length(stat)]){
                        v <- i+fDist[y,6]
                        u <- i+(fDist[y,6]+1)
                        sess <- i+fDist[y,6]
                        z[i] <- Dataset$FData[y]*r
                } else {
                    if(stat[2] > stat[length(stat)])
                        v <- i+(fDist[y,6]+1)
                    u <- i+fDist[y,6]
                    sess <- i+fDist[y,6]
                    z[i] <- Dataset$FData[y]*r
        }
        r <-((1-exp(-(u/v)^k))/(1-exp(-(v/v)^k)))
        z[i] <- z[i-1]*r
        Pred[i,] = c(z[i])
}

                Pr <- format(round(Pred, 2), nsmall = 2)
                cat("Predictions, significant digits = 2, rows 1:20:","\n", "\n")
                                print(Pr[1:20])

                dup <- format(round(Pred, 2), nsmall = 2)
                (xu1 <- dup[duplicated(dup, fromLast = TRUE)])
                        xu2 <- xu1[1]
cat("\n")
cat("First predicted duplicate value")
print(xu2[1])
dup2 <- which(dup[,1] == xu2[1])
saveRDS(dup2, "dup2.rds")
cat("\n")
cat("Row numbers for predicted sessions")
print(dup2)
cat("\n")
cat("Probability curve for SPRE predictions, rows 1:6:")
cat("\n")

#' @name x
#' @description The generic matrix function creates the first duplicated point estimations. The associated session is then the effective clinical stability.
#' @param m2 A number
#' @param pv A number
#' @return plot
#' @export
#'
#'  Computation for effective clinical stability using point estimation of Weibull CDF

                x <- Pred
                m2 <- readRDS("qn.rds")
                pv <- readRDS("pv.rds")
                pv <- (signif(pv,digits=2))

        par(mfrow = c(1,1), mar=c(5.5,5.2,5.5,3.0))
                plot(x, main = "Probability for SPRE Predictions \n with
     effective clinical stability(ECS)",
                xlab = (paste("Prediction & Pvalue from Change Point=",m2,", p=", pv)),
                ylab = "Probability",
                xlim = c(0,40), ylim = c(0,1))
    
cat("\n")
cat("Clinical stability is determined from the probability curve for SPRE predictions,if any 2 predictions from the second significant digit (after the decimal place) remains constant. Otherwise, there is clinical significance at the vertical line (zero) at the change point.")

#'  Computation for 'wiebull' cumulative probability distribution

        pweibull(x, shape=k, scale = sL, lower.tail = TRUE, log.p = FALSE)
                curve(pweibull(x, scale=sL, shape=k),
                from=2, to=40, add=TRUE, lwd=2)

cat("\n","\n")
cat("First six predicted probability values for Weibull distrbution","\n")

        m4 <- print(pweibull(x, shape=k, scale = sL, lower.tail = TRUE, log.p               =FALSE)[1:6])
                b <- abline(v=0, col="blue", lty=3)
                g <-if (dup2[[2]] <- 0) then (dup2[[2]] <- 1:1)
                        rm(b,g)

                m3 <- readRDS("dup2.rds")
                c <- 1
                a <- abline(v=m3, col="red", lty=2, lwd=1.5)
                y <- if (c==1) a else b

cat("\n", "\n")

                dev.set(dev.next())
                rm (list = ls())
               
