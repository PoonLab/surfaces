setwd("~/git/surfaces/data")

dnds <- read.csv("step7_dnds.csv")
mdat <- read.csv("metadata.csv")

idx <- match(paste(dnds$virus, dnds$protein), 
             paste(mdat$virus, mdat$protein)) 

dnds <- cbind(dnds, mdat[idx, !is.element(names(mdat), c('virus', 'protein'))])

# alphabetical order
x <- as.integer(as.factor(paste(dnds$family, dnds$virus)))

pdf("dnds.pdf", width=6, height=5)
par(mar=c(5,5,5.5,1))
plot(x, dnds$dnds, log='y', xaxt='n', bty='n', las=2, tck=-0.02,
     xlab=NA, ylab="mean dN / dS", cex.axis=0.8, mgp=c(3,0.5,0),
     cex=ifelse(dnds$enveloped, 1, 0.9),
     pch=ifelse(dnds$enveloped, 
                ifelse(dnds$exposed, 19, 1), 
                ifelse(dnds$exposed, 23, 5)),
     col=ifelse(dnds$enveloped, 'darkorange2', 'cadetblue'),
     bg='cadetblue')
labels <- unique(dnds$abbrv[order(x)])
axis(side=1, at=1:length(labels), labels=labels, las=2, cex.axis=0.8)

counts <- sapply(split(dnds$virus, dnds$family), function(x) length(unique(x)))
right <- cumsum(counts)
left <- right-counts
axis(side=3, at=left + (right-left+1)/2, labels=names(fams), 
     las=2, cex.axis=0.7, line=-0.5, lwd=0)
is.env <- sapply(split(dnds$enveloped, dnds$family), function(x) x[1])
segments(x0=left+0.75, x1=right+0.25, y0=1, xpd=NA, 
         col=ifelse(is.env, 'darkorange2', 'cadetblue'), lwd=2)
#legend(x=0.5, y=0.005, legend=c("surface-exposed"), pch=19, bty='n',
#       x.intersp=0.8, cex=0.8, pt.cex=1.1, y.intersp=0, col='grey50',
#       text.col='grey50')
dev.off()


hist(log(dnds$dnds))

# naive model
fit <- lm(log(dnds$dnds) ~ dnds$exposed)
summary(fit)
confint(fit)

fit <- glm(dnds ~ exposed, data=dnds, family=Gamma(link='log'))
summary(fit)
confint(fit)

require(lme4)
fit <- lmer(log(dnds$dnds) ~ dnds$exposed + dnds$enveloped + (1|dnds$virus))
summary(fit)
confint(fit)

require(BSDA)
EDA(residuals(fit))


fit2 <- glmer(dnds ~ exposed  + (1|virus), data=dnds, family=Gamma(link='log'))
summary(fit2)
EDA(residuals(fit2))
confint(fit2, method='boot', boot.type='basic')
