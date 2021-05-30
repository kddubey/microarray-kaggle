
# for internal use
# Install the penalizedSVM package into my Kaggle workspace
install.packages("penalizedSVM", lib="/kaggle/working")
library(penalizedSVM, lib="/kaggle/working")

# Modeling
library(penalizedSVM) # scadsvc fxn
library(caret)        # confusionMatrix fxn
library(plyr)         # ldply fxn

# Cleaning
library(tidyverse)    # select (masks MASS::select), spread fxns
library(reshape2)     # melt fxn

PATH.TO.DATA = '../input/gene-expression'

PATH.TO.TRAIN = file.path(PATH.TO.DATA, 
                          'data_set_ALL_AML_train.csv')
PATH.TO.TEST  = file.path(PATH.TO.DATA, 
                          'data_set_ALL_AML_independent.csv')

PATH.TO.LABELS = file.path(PATH.TO.DATA, 'actual.csv')

# X.tr will represent the feature matrix for the training set
X.tr = read.csv(PATH.TO.TRAIN,
                stringsAsFactors=FALSE, quote="")
# X.te will represent the feature matrix for the test set
X.te = read.csv(PATH.TO.TEST,
                stringsAsFactors=FALSE, quote="")

y = read.csv(PATH.TO.LABELS,
             stringsAsFactors=FALSE, quote="")

head(X.tr)[,1:8]

print("Patients in the training set:")
colnames(dplyr::select(X.tr, contains("X")))

print("Patients in test set:")
colnames(dplyr::select(X.te, contains("X")))

print("Patients in labels dataset are in order:")
print(!is.unsorted(y[,"patient"]))

clean.gene.data = function(X) {
  # This code is nearly identical to Nelson Gonzabato's notebook:
  #   kaggle.com/gonnel/who-is-at-risk-of-cancer-a-simple-analysis
  # I've only added explanation and put it in a function among
  # other minor changes. 
  
  # I'm going to ignore the "call" column. 
  # See the discussion here to learn about it:
  #   kaggle.com/crawford/gene-expression/discussion/120087
  # The "Gene Description" column is redundant for our
  # purposes so we'll remove that column as well. 
  X = dplyr::select(X, -contains("call"), -Gene.Description)
  
  # Sort examples by patient number so that they are in the same
  # order as the labels (they are originally in a different
  # order than the labels!):
  X = X[, str_sort(names(X), numeric=TRUE)]
  
  # Convert X into a standard feature matrix so that X[i,j] is
  # the value of gene expression j for patient i. 
  X = melt(X, id.vars="Gene.Accession.Number")
  X = spread(X, Gene.Accession.Number, value)
  X = dplyr::select(X, -variable)
  
  # Note that features were already re-scaled, which is necessary
  # for SVMs. 
  
  return(X)
}

X.tr = clean.gene.data(X.tr)
X.te = clean.gene.data(X.te)

ncol(X.tr)

X.tr = as.matrix(X.tr)
X.te = as.matrix(X.te)

y.vec = y[,"cancer"]
y.bin = rep(0, length(y.vec))
y.bin[y.vec == "ALL"] = 1
y.bin[y.vec == "AML"] = -1

n.tr = nrow(X.tr) 
n.te = nrow(X.te)

n = n.tr + n.te

y.tr = y.bin[1:n.tr]
y.te = y.bin[(n.tr+1):n]

# You can also write these data for future use so that you don't
# have to clean every time. 
#write.csv(X.tr, "train/features.csv", row.names=FALSE)
#write.csv(X.te, "test/features.csv",  row.names=FALSE)
#write.csv(data.frame(label=y.tr), "train/labels.csv",
#          row.names=FALSE)
#write.csv(data.frame(label=y.te), "test/labels.csv",
#          row.names=FALSE)

# Then, in a new R session you'd run these lines:
#X.tr = as.matrix(read.csv("golub/train/features.csv"))
#y.tr = read.csv("golub/train/labels.csv")$label
#X.te = as.matrix(read.csv("golub/test/features.csv"))
#y.te = read.csv("golub/test/labels.csv")$label

print(paste(ncol(X.tr), "features (genes),",
            nrow(X.tr), "examples (patients)"))

tune.hp.loocv = function(X, y, lambdas, iters, verbose=TRUE) {
  # Returns a matrix called errors where errors[i,j] is the
  # proportion of LOOCV model misclassifications for the SCAD SVM
  # trained using hyperparameter settings lambda=lambdas[i] and
  # maxIter=iters[j].
  # If verbose, then the number of misclassifications is printed
  # along with the lambda and maxIter setting. By default
  # verbose=TRUE. 
  
  # Implementation notes:
  # - The try clause handles the case where no genes are selected.
  #   This occurs when lambda and maxIter are set too high. 
  # - The seed of the SCAD SVM optimizer is 123 by default. 
  
  m = length(lambdas)
  p = length(iters)
  errors = matrix(rep(1, m*p), m, p)
  n = length(y)
  
  lowest.err = n
  best.svm.outputs = rep(NA, n)
  
  for (i in 1:m) {
    lambda = lambdas[i]
    for (j in 1:p) {
      iter = iters[j]
      err = 0
      num.models = 0 # number of models where a gene is selected
      svm.outputs = rep(NA, n)
      for (k in 1:n) {
        X.tr.fold = X[-k,]
        x.val = X[k,]
        y.tr.fold = y[-k]
        y.val = y[k]
        model = scadsvc(X.tr.fold, y=y.tr.fold, lambda=lambda,
                        maxIter=iter, verbose=FALSE)
        try({
          svm.out = x.val[model$xind] %*% model$w + model$b
          svm.outputs[k] = svm.out
          pred = sign(svm.out)
          err = err + (pred != y.val)
          num.models = num.models + 1
        })
      }
      if (verbose) {
        print(paste("lambda =", lambda, "iter =", iter,
                    "misclassified", err, "points"))
      }
      if (err < lowest.err) {
        lowest.err = err
        best.svm.outputs = svm.outputs
      }
      errors[i,j] = err / num.models
    }
  }
  return(list("errors"=errors, "outputs"=best.svm.outputs))
}

# Run and time HP search (takes ~12 minutes in the Kaggle engine)
lambdas = seq(0, 1.5, by=0.25) # random search: runif(7, 0, 1.5)
iters = c(1000, 2000)

tune.time = system.time({
  loocv.result = tune.hp.loocv(X.tr, y.tr, lambdas, iters)
})

errors = loocv.result$errors

print(paste("Tuning took", round(tune.time[["elapsed"]]/60, 2),
            "minutes."))

plot.hp.error = function(hp1, hp2, errors, hp1.name, hp2.name) {
    # Plots error according to hyperparemter settings. 
    # hp1 appears on the x-axis. hp2 is paritioned into
    # different colors. 
    # errors is the matrix such that hyperparameter settings
    # hp1[i] and hp2[j] resulted in errors[i,j]. 
    plot(hp1, errors[,1], type="l",
         xlab=hp1.name, ylab="Proportion misclassified",
         main="LOOCV error of classifier",
         ylim=c(0, max(errors)+0.1))
    for (i in 2:length(hp2)) {
        lines(hp1, errors[,i], col=i)
    }
    legend("topleft", legend=paste(hp2.name, "=", hp2),
           col=1:length(hp2), box.lty=0, lty=c(1,1), cex=0.75)
}

plot.hp.error(lambdas, iters, errors, 
              "lambda", "maxIters")

plot.hp.error(iters, lambdas, t(errors), 
              "maxIters", "lambda")

# Here's an interactive 3D plot that you can run for "fun"
#library(plotly) # plot_ly, layout functions
#errors.df = data.frame(x = rep(lambdas, each=length(iters)),
#                       y = rep(iters, length(lambdas)),
#                       z = as.vector(t(errors)))
#axx = list(title = "lambda")
#axy = list(title = "maxIters")
#axz = list(title = "Number misclassified")
#hp.fig = plot_ly(x=~errors.df$x, y=~errors.df$y,
#                 z=~errors.df$z, type='mesh3d')
#hp.fig = hp.fig %>% layout(scene = list(xaxis=axx, yaxis=axy, 
#                                        zaxis=axz))
#hp.fig

best.ind = which(errors == min(errors), arr.ind=TRUE)
lambda = lambdas[best.ind[[1]]] # 1.25
maxIter = iters[best.ind[[2]]]  # 2000

# The seed of the scadsvc function is set to 123 by default
model = scadsvc(X.tr, y=y.tr, lambda=lambda, maxIter=maxIter,
                verbose=TRUE)

preds.tr = sign(X.tr[,model$xind] %*% model$w + model$b)
sum(preds.tr == y.tr)

preds.te = sign(X.te[,model$xind] %*% model$w + model$b)
sum(preds.te == y.te)

# which patient in the test set was misclassified?
misclassifieds = preds.te != y.te
paste("Patient", which(misclassifieds) + n.tr, 
      "was misclassified as label", preds.te[misclassifieds])

colnames(X.tr)[model$xind]

bootstrap.data = function(X, y, 
                          N.sim=20, lambda=1.25, maxIter=2000) {
    # model.coeffs[i,j] is the coefficient estimate for gene i
    # in simulation j
    n = nrow(X)
    d = ncol(X)
    model.coeffs = matrix(rep(0, d*N.sim), d, N.sim)
    # count.selected[i] is the number of times that gene i is
    # selected over N.sim simulations
    count.selected = rep(0, d)
    for (k in 1:N.sim) {
      # form bootstrap sample
      set.seed(k*123)
      boot.inds = sample(n, n, replace=TRUE)
      X.boot = X[boot.inds,]
      y.boot = y[boot.inds]
      # train model on this data
      model.boot = scadsvc(X.boot, y=y.boot, lambda=lambda,
                           maxIter=maxIter, verbose=FALSE)
      if (model.boot == 'No variable selected.') next
      # update variables
      model.coeffs[model.boot$xind, k] = model.boot$w
      for (ind in model.boot$xind) {
          count.selected[ind] = count.selected[ind] + 1
      }
    }
    return(list("model.coeffs"=model.coeffs,
                "count.selected"=count.selected))
}

plot.selected = function(count.selected, N.sim, 
                         count.thresh=5) {
    # Plot number of simulations that each gene was selected,
    # i.e., has non-zero weight. 
    # Only plot genes that are selected at least count.thresh
    # times.  
    selected.inds = which(count.selected >= count.thresh)
    selected.genes = colnames(X.tr)[selected.inds]
    num.selections = count.selected[selected.inds]
    # plot
    title = paste("Genes selected", count.thresh, 
                  "times or more in", 
                  N.sim, "bootstrapped datasets")
    ggplot(mapping=aes(x=reorder(selected.genes,-num.selections), 
                       y=num.selections)) +
        geom_bar(stat="identity", fill="steelblue") +
        theme(axis.text.x=element_text(angle=90),
              plot.title=element_text(face='bold',hjust=0.5)) +
        labs(title=title,
             x="Gene accession number", 
             y="Number of simulations selected")
}

X = rbind(X.tr, X.te)
y = c(y.tr, y.te)

N.sim = 100 # number of bootstrap simulations
lambda = 1.25
maxIter = 2000

bootstrap.result = bootstrap.data(X, y, N.sim, 
                                  lambda=lambda, maxIter=maxIter)

model.coeffs = bootstrap.result$model.coeffs
count.selected = bootstrap.result$count.selected
plot.selected(count.selected, N.sim)

sim.counter = table(colSums(model.coeffs != 0))
sim.df = as.data.frame(sim.counter)
colnames(sim.df) = c("# genes selected", "# simulations")
sim.df

sum(count.selected != 0)

top4.genes = c('M19507_at', 'M27891_at', 'M96326_rna1_at', 
               'Y00787_s_at')
top4.gene.inds = match(top4.genes, colnames(X))
num.together = colSums(model.coeffs[top4.gene.inds,] != 0)
as.data.frame(table(num.together))

# form dataframe
top4.df = data.frame(t(model.coeffs[top4.gene.inds,]))
names(top4.df) = top4.genes
top4.df = gather(top4.df, "gene", "coefficient")
# only examine coefficient if it's non-zero
top4.df = top4.df[top4.df$coefficient != 0, ]

for (gene in top4.genes) {
    gene.df = top4.df[top4.df$gene == gene,]
    p = ggplot(gene.df, aes(x=coefficient, 
                            color=gene, 
                            fill=gene)) +
            geom_histogram(alpha=1/3, bins=5)
    print(p)
}

mean(top4.df$coefficient < 0)

min(top4.df$coefficient)

max(top4.df$coefficient)

print(model$w)

calibration.stats = function(labels, pred.probs, bins=10) {
    # bins is an upper bound on the number of bins used to
    # measure calibration
    require(plyr)
    bin.pred = cut(pred.probs, bins)
    res = ldply(levels(bin.pred), function(x) {
        idx = (x == bin.pred)
        pred.mean = mean(pred.probs[idx])
        obs.mean = mean(labels[idx])
        bin.size = sum(idx)
        se = sqrt((obs.mean * (1 - obs.mean)) / bin.size)
        ll = max(0, obs.mean - 1.96*se)
        ul = min(1, obs.mean + 1.96*se)
        c(pred.mean, obs.mean, bin.size, ll, ul)
    })
    # remove bins with no data
    colnames(res) = c("pred.mean", "obs.mean", "bin.size", 
                      "ll", "ul")
    is.nan.idx = !is.nan(res$pred.mean)
    res = res[is.nan.idx,]
    return(res)
}

plot.calibration = function(calibration.res, 
                            title="Calibration curve") {
    ggplot(res, aes(x=pred.mean, y=obs.mean, 
                    color="Model calibration")) +
    geom_abline(aes(slope=1, intercept=0, 
                    color="Perfect calibration")) +
    geom_point() +
    geom_text(aes(label=bin.size), hjust=2, vjust=0) +
    geom_errorbar(aes(ymin=ll, ymax=ul), width=0) +
    theme(plot.title = element_text(face="bold", hjust=0.5)) +
    labs(title = title, 
         x = "Mean predicted probability of ALL", 
         y = "Observed fraction of ALL",
         color="")
}

outputs.loocv = loocv.result$outputs
sum(is.na(outputs.loocv))

outputs.tr = X.tr[,model$xind] %*% model$w + model$b
na.inds = which(is.na(outputs.loocv))
outputs.loocv[na.inds] = outputs.tr[na.inds]

transform.target = function(y) {
  # Transforms the binary (0,1 or -1,1) target into
  # MAP-estimated, soft targets as described in Platt (1999)
  y.new = rep(NA, length(y))
  is.pos = (y == 1)
  is.neg = (y != 1)
  n.pos = sum(is.pos)
  n.neg = sum(is.neg)
  y.new[which(is.pos)] = (n.pos + 1) / (n.pos + 2)
  y.new[which(is.neg)] = 1 / (n.neg + 2)
  return(y.new)
}

data.tr = data.frame(x=outputs.loocv, y=transform.target(y.tr))

# optimization functions specifically for Platt scaling
pred.prob = function(x, par) {
  X = cbind(rep(1, length(x)), x)
  z = X %*% par
  return(1 / (1 + exp(-z)))
}

kld.loss = function(par, data) {
  x = data$x
  y = data$y
  p = pred.prob(x, par)
  return(-sum(y*log(p) + (1-y)*log(1-p)))
}

kld.loss.grad = function(par, data) {
  x = data$x
  y = data$y
  p = pred.prob(x, par)
  X = cbind(rep(1, length(x)), x)
  return(-t(X) %*% (y-p))
}

# perform optimization
# no need for L-BFGS as this is just a 2-parameter problem
init.par = c(0, 0) # intial intercept, slope
result = optim(init.par, kld.loss,
               gr=kld.loss.grad, data=data.tr,
               method='BFGS')

outputs.te = X.te[,model$xind] %*% model$w + model$b
pred.te.probs = pred.prob(outputs.te, result$par)

hist(pred.te.probs, breaks=15, 
     xlab="Predicted probability of ALL",
     ylab="Frequency",
     main="Predicted probabilities in test set")

title = "Calibration after regularization"
y.te.sparse = pmax(0, y.te)
res = calibration.stats(y.te.sparse, pred.te.probs, bins=5)
plot.calibration(res, title=title)

data.tr = data.frame(x=outputs.tr, y=pmax(0, y.tr))
model.calibrated = glm(y ~ x, family=binomial, data=data.tr)

pred.te.probs.base = predict(model.calibrated, 
                             data.frame(x=outputs.te),
                             type='response')
hist(pred.te.probs.base, breaks=15, 
     xlab="Predicted probability of ALL",
     ylab="Frequency",
     main="Predicted probabilities in test set")

title = "Calibration without regularization"
res = calibration.stats(y.te.sparse, pred.te.probs.base, bins=5)
plot.calibration(res, title=title)

pred.te.probs[misclassifieds]

pred.te.probs.base[misclassifieds]

conf.mat = confusionMatrix(as.factor(preds.te), 
                           as.factor(y.te),
                           positive="1")
conf.mat$table

tn = conf.mat$table[1,1] # true negatives
fn = conf.mat$table[1,2] # false negatives
fp = conf.mat$table[2,1] # false positives
tp = conf.mat$table[2,2] # true positives

# Pr(patient has ALL | model predicts patient has ALL) ≈
ppv = tp/(tp + fp)
# Pr(patient has AML | model predicts patient has AML) ≈
npv = tn/(tn + fn)

# Pr(model predicts patient has ALL | patient has ALL) ≈
tpr = tp/(tp + fn)
# Pr(model predicts patient has AML | patient has AML) ≈
tnr = tn/(tn + fp)

metrics = c("Positive predictive value (PPV)", 
            "Negative predictive value (NPV)", 
            "True positive rate (TPR)", 
            "True negative rate (TNR)")
data.frame(metric=metrics, score=c(ppv, npv, tpr, tnr))

hoeffding.inequality = function(m, n, eps) {
    return( m * 2*exp(-2*n*eps^2) )
}

m = 2
eps = 0.1
hoeffding.inequality(m, n.te, eps)

hoeffding.sample.size = function(m, p, eps) {
    (log(2) + log(m) - log(p)) / (2*eps^2)
}

m = 1
eps = 0.1
p = 0.1
hoeffding.sample.size(m, p, eps)
