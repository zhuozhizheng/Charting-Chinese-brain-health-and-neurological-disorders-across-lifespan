library(gamlss)

fit_model=function(num){
  
  i_poly<-list_par[num,1]
  j_poly<-list_par[num,2]
  i_rnd<-list_par[num,3]
  j_rnd<-list_par[num,4]
  
  list_fit<-data.frame(matrix(0,1,9));
  colnames(list_fit)<-c('mu_poly','sigma_poly','mu_random','sigma_random','Global_Deviance','AIC','SBC','BIC','Selection');
  
  if(i_rnd==0&j_rnd==0)
  {m0<-gamlss(formula=feature~fp(Age,npoly=i_poly)+Sex,
              sigma.formula =feature~fp(Age,npoly=j_poly)+Sex,
              #nu.formula = feature~1,
              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),control=con,
              data=data1)

  } else if(i_rnd==0&j_rnd==1)
  {m0<-gamlss(formula=feature~fp(Age,npoly=i_poly)+Sex,
              sigma.formula =feature~fp(Age,npoly=j_poly)+Sex+random(Site_ZZZ),
              #nu.formula = feature~1,
              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),control=con,
              data=data1)
  
  
  } else if(i_rnd==1&j_rnd==0){m0<-gamlss(formula=feature~fp(Age,npoly=i_poly)+Sex+random(Site_ZZZ),
                                          sigma.formula =feature~fp(Age,npoly=j_poly)+Sex,
                                          #nu.formula = feature~1,
                                          family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),control=con,
                                          data=data1)
  
  
  } else if(i_rnd==1&j_rnd==1){m0<-gamlss(formula=feature~fp(Age,npoly=i_poly)+Sex+random(Site_ZZZ),
                                          sigma.formula =feature~fp(Age,npoly=j_poly)+Sex+random(Site_ZZZ),
                                          #nu.formula = feature~1,
                                          family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),control=con,
                                          data=data1)
  
  
  } 
  m1<-m0
  

  num=1;
  list_fit[num,'mu_poly']<-i_poly;list_fit[num,'sigma_poly']<-j_poly;list_fit[num,'mu_random']<-i_rnd;list_fit[num,'sigma_random']<-j_rnd;
  
  list_fit[num,'Global_Deviance']<-m1$G.deviance;list_fit[num,'AIC']<-m1$aic;list_fit[num,'SBC']<-m1$sbc;list_fit[num,'BIC']<-BIC(m1);
  list_fit[num,'Selection']<-num;
  
  list_fit1<-list_fit;
  return(list_fit1)
  
}



best_fit=function(sel_mu_poly,sel_sigma_poly,i_rnd,j_rnd){
  
  
  
  if(i_rnd==0&j_rnd==0)
  {m0<-gamlss(formula=feature~fp(Age,npoly=sel_mu_poly)+Sex,
              sigma.formula =feature~fp(Age,npoly=sel_sigma_poly)+Sex,
              #nu.formula = feature~1,
              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),control=con,
              data=data1)

 
  } else if(i_rnd==0&j_rnd==1)
  {m0<-gamlss(formula=feature~fp(Age,npoly=sel_mu_poly)+Sex,
              sigma.formula =feature~fp(Age,npoly=sel_sigma_poly)+Sex+random(Site_ZZZ),
              #nu.formula = feature~1,
              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),control=con,
              data=data1)
  

  
  } else if(i_rnd==1&j_rnd==0){m0<-gamlss(formula=feature~fp(Age,npoly=sel_mu_poly)+Sex+random(Site_ZZZ),
                                          sigma.formula =feature~fp(Age,npoly=sel_sigma_poly)+Sex,
                                          #nu.formula = feature~1,
                                          family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),control=con,
                                          data=data1)
  
  
  } else if(i_rnd==1&j_rnd==1){m0<-gamlss(formula=feature~fp(Age,npoly=sel_mu_poly)+Sex+random(Site_ZZZ),
                                          sigma.formula =feature~fp(Age,npoly=sel_sigma_poly)+Sex+random(Site_ZZZ),
                                          #nu.formula = feature~1,
                                          family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),control=con,
                                          data=data1)
  
  
  } 
  
  return(m0)
  
}

zzz_cent<-function(obj=m2,type=c("centiles","z-scores"),mu=res$mu,sigma=res$sigma,nu=res$nu,
                   cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age,yval=NULL,
                   calibration=FALSE,lpar=3,cdf=NULL)
{
  calc.cent <- function(xvar, cent,mu,sigma,nu,lpar) {
    o <- order(xvar)
    mat <- xvar[o]
    cent <- cent
    for (var in cent) {
      if (lpar == 1) {
        newcall <- call(qfun, var/100, mu = mu[o])
      }
      else if (lpar == 2) {
        newcall <- call(qfun, var/100, mu = mu[o], sigma = sigma[o])
      }
      else if (lpar == 3) {
        newcall <- call(qfun, var/100, mu = mu[o], sigma = sigma[o], 
                        nu = nu[o])
      }
      else {
        newcall <- call(qfun, var/100, mu = mu[o], sigma = sigma[o], 
                        nu = nu[o], tau = tau[o])
      }
      ll <- eval(newcall)
      mat <- cbind(mat, ll)
    }
    mat <- as.data.frame(mat)
    nnn <- paste("C", as.character(cent), sep = "")
    names(mat) <- c(xname, nnn)
    return(mat)
  }
  
  
  if (type == "centiles") {
    fname <- obj$family[1]
    qfun <- paste("q", fname, sep = "")
    xvar <- xvalues
    if (calibration) {
      z <- quantile(resid(obj), probs = cent/100)
      p <- pNO(z, mu = 0, sigma = 1)
      cent <- round(100 * p, digits = 2)
    }
    mat <- calc.cent(xvar = xvar, cent = cent,lpar=lpar,mu=mu,sigma=sigma,nu=nu)
    colnames(mat) <- c("x", as.character(cent))
    p2<-mat;
    return(p2)
  }
  
  if (type == "z-scores"&is.null(cdf)) {
    if (calibration)
      stop("calibration is not implemeted yet in z-scores")
    if (is.null(yval))
      stop("the y values should be set if type=z-scores is used")
    if (length(yval) != length(xvalues))
      stop("length of xvalues and yval is not the same")
    fname <- obj$family[1]
    qfun <- paste("p", fname, sep = "")
    if (lpar == 1) {
      newcall <- call(qfun, yval, mu = mu)
    }
    else if (lpar == 2) {
      newcall <- call(qfun, yval, mu = mu, sigma = sigma)
    }
    else if (lpar == 3) {
      newcall <- call(qfun, yval, mu = mu, sigma = sigma,
                      nu = nu)
    }
    else {
      newcall <- call(qfun, yval, mu = mu, sigma = sigma,
                      nu = nu, tau = tau)
    }
    cdf <- eval(newcall)
    rqres <- qnorm(cdf)
    return(rqres)
  } else if (type == "z-scores"&!is.null(cdf)) {
    if (calibration)
      stop("calibration is not implemeted yet in z-scores")
    if (is.null(yval))
      stop("the y values should be set if type=z-scores is used")
    if (length(yval) != length(xvalues))
      stop("length of xvalues and yval is not the same")
    fname <- obj$family[1]
    qfun <- paste("p", fname, sep = "")
    if (lpar == 1) {
      newcall <- call(qfun, yval, mu = mu)
    }
    else if (lpar == 2) {
      newcall <- call(qfun, yval, mu = mu, sigma = sigma)
    }
    else if (lpar == 3) {
      newcall <- call(qfun, yval, mu = mu, sigma = sigma,
                      nu = nu)
    }
    else {
      newcall <- call(qfun, yval, mu = mu, sigma = sigma,
                      nu = nu, tau = tau)
    }
    cdf <- eval(newcall)
    rqres <- cdf
    return(rqres)
  }
  

}



