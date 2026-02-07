#' @title Operating Characteristics for the Inserted Shared Keyboard Design
#' 
#' @description
#' A short description...
#' 
#' @param target_prob description
#' 
#' @export

get_OC_Insert_SKBD <- function(
    target_prob, tox_prob, 
    n_cohort, cohort_size,
    dose_set, dose_range,
    C1 = 0.6, C2 = 0.6, 
    start_dose = 1,
    margin_left = 0.05, margin_right = 0.05, 
    n_earlystop = 1000, 
    cutoff_elimin = 0.95,
    extra_safe = FALSE, offset = 0.05, 
    symmetric = FALSE, theta = NULL,
    shared = TRUE, # incorporate the keybooard design as a special case
    light_return = TRUE,
    n_trial = 1000, seed = 6
) {
  set.seed(seed)
  
  ## ---- argument checks ----
  # target
  if (!is.numeric(target_prob) || length(target_prob) != 1 || is.na(target_prob) ||
      target_prob <= 0 || target_prob >= 1) {
    stop("`target_prob` must be a single number in (0, 1).")
  }
  if (target_prob < 0.05) warning("`target_prob` is very low (< 0.05); OC may be unstable.")
  if (target_prob > 0.5)  warning("`target_prob` is very high (> 0.5); uncommon in phase I.")
  
  # margins
  if (!is.numeric(margin_left) || !is.numeric(margin_right) ||
      length(margin_left) != 1 || length(margin_right) != 1 ||
      anyNA(c(margin_left, margin_right)) || margin_left < 0 || margin_right < 0) {
    stop("`margin_left` and `margin_right` must be nonnegative scalars.")
  }
  key_L <- target_prob - margin_left
  key_U <- target_prob + margin_right
  if (key_L <= 0 || key_U >= 1 || key_L >= key_U) {
    stop("Invalid target key: require 0 < target_prob - margin_left < target_prob + margin_right < 1.")
  }
  
  # tox_prob
  if (!is.numeric(tox_prob) || anyNA(tox_prob) || any(tox_prob < 0 | tox_prob > 1)) {
    stop("`tox_prob` must be numeric in [0, 1] with no NA.")
  }
  
  # n_cohort / cohort_size / n_trial
  int1 <- function(x) is.numeric(x) && length(x) == 1 && !is.na(x) && x %% 1 == 0 && x >= 1
  if (!int1(n_cohort))    stop("`n_cohort` must be a positive integer.")
  if (!int1(cohort_size)) stop("`cohort_size` must be a positive integer.")
  if (!int1(n_trial))     stop("`n_trial` must be a positive integer.")
  
  # n_earlystop (here used as posterior-sample size option placeholder; we keep the check)
  if (!int1(n_earlystop)) stop("`n_earlystop` must be a positive integer.")
  if (n_earlystop <= 6) {
    warning("`n_earlystop` is small; Recommended range is 9-18.")
  }
  
  # cutoff / extra safety
  if (!is.numeric(cutoff_elimin) || length(cutoff_elimin) != 1 || is.na(cutoff_elimin) ||
      cutoff_elimin <= 0 || cutoff_elimin >= 1) {
    stop("`cutoff_elimin` must be a scalar in (0, 1).")
  }
  if (!is.numeric(offset) || length(offset) != 1 || is.na(offset) || offset < 0 || offset >= 0.5) {
    stop("`offset` must be a scalar in [0, 0.5).")
  }
  if (extra_safe && (cutoff_elimin - offset) <= 0) {
    stop("When `extra_safe=TRUE`, require `cutoff_elimin - offset > 0`.")
  }
  
  # dose_set vs tox_prob
  if (!is.numeric(dose_set) || anyNA(dose_set)) stop("`dose_set` must be numeric with no NA.")
  if (length(dose_set) != length(tox_prob)) {
    stop("`dose_set` must have the same length as `tox_prob` (i.e., n_dose).")
  }
  if (anyDuplicated(dose_set)) stop("`dose_set` must have unique values (no duplicates).")
  if (is.unsorted(dose_set, strictly = TRUE)) {
    stop("`dose_set` must be strictly increasing (sorted).")
  }
  
  # dose_range interval
  if (!is.numeric(dose_range) || length(dose_range) != 2 || anyNA(dose_range)) {
    stop("`dose_range` must be a numeric vector of length 2: c(lower, upper).")
  }
  if (dose_range[1] >= dose_range[2]) stop("`dose_range[1]` must be < `dose_range[2]`.")
  if (min(dose_set) < dose_range[1] - 1e-12 || max(dose_set) > dose_range[2] + 1e-12) {
    stop("All `dose_set` values must lie inside `dose_range`.")
  }
  
  # start_dose index (in prespecified doses)
  if (!int1(start_dose) || start_dose < 1 || start_dose > length(tox_prob)) {
    stop("`start_dose` must be a single integer in 1:n_dose (indexing prespecified doses).")
  }
  
  # theta / kernel mode
  # Standardize dose to [0,1] using dose_range (stabilizes theta meaning)
  dose_set_std = (dose_set - dose_range[1]) / (dose_range[2] - dose_range[1])
  if (is.null(theta)) {
    dose_diff_min = min(diff(dose_set_std))
    # defaults (safe-ish): asym borrows more from higher doses => theta1 > theta2
    theta <- if (symmetric) {
      -log(0.5) / dose_diff_min^2
    } else {
      c(-log(0.2), -log(0.8)) / dose_diff_min^2
    }
  }
  if (symmetric) {
    if (!is.numeric(theta) || length(theta) != 1 || is.na(theta) || theta <= 0) {
      stop("When `symmetric=TRUE`, `theta` must be a single positive number.")
    }
  } else {
    if (!is.numeric(theta) || length(theta) != 2 || anyNA(theta) || any(theta <= 0)) {
      stop("When `symmetric=FALSE`, `theta` must be a length-2 positive vector c(theta1, theta2).")
    }
    if (!(theta[1] > theta[2])) {
      stop("For asymmetric kernel, require theta[1] > theta[2] (borrow less from lower doses).")
    }
  }
  
  
  n_dose     = length(tox_prob)       # number of dose levels
  n_patients = n_cohort * cohort_size # total number of patients
  
  sel_dose_idx     = integer(n_trial) # selected dose level index (j*)
  sel_dose         = numeric(n_trial) # selected dose value (d_j*)
  insert_at_cohort = integer(n_trial) # cohort index when first insertion occurs (0 if none)
  n_insertions     = integer(n_trial) # number of insertions in each trial

  pri_alpha = pri_beta = 0.5 # hyperparameters for beta prior
  sim_data = NULL           # store simulated output (will rbind/append)
  
  ## get keys
  key_cutpoints  = get_Key(target_prob, margin_left, margin_right)
  target_key_idx = which.min(abs(key_cutpoints - (target_prob - margin_left)))
  
  
  ##----------------Simulate trials-----------------##
  for(trial in 1:n_trial){
    d = start_dose  # current dose index (starting dose for this trial)
    
    # --- working dose grid (may expand after insertion) ---
    dose_set_work = dose_set_std   # current dose set; starts with prespecified doses and may add inserted doses
    tox_prob_work = tox_prob       # true toxicity probs aligned with dose_set_work; expanded when a new dose is inserted
    n_dose_work   = length(tox_prob_work)  # current number of dose levels on the working grid
    
    # --- trial-level sufficient statistics on the working grid ---
    n_dlt_work     = rep(0, n_dose_work)   # number of DLTs observed at each dose level
    n_treated_work = rep(0, n_dose_work)   # number of patients treated at each dose level
    
    # --- BKP posterior parameters (initialized at Beta prior) ---
    post_alpha_work = rep(pri_alpha, n_dose_work)  # posterior alpha parameters at each dose level
    post_beta_work  = rep(pri_beta,  n_dose_work)  # posterior beta  parameters at each dose level
    
    # --- safety tracking ---
    is_eliminated = rep(FALSE, n_dose_work)  # dose elimination indicator (TRUE = eliminated / not allowed)
    is_earlystop  = FALSE                   # early termination flag (TRUE = stop the trial for safety)
    
    # --- insertion bookkeeping ---
    n_insertions_trial = 0  # number of dose insertions performed in this trial
    
    
    # --- main trial loop: treat cohorts until reaching total sample size or early stop ---
    while(TRUE){
      # if enough patients have been treated at dose d, exit the loop
      if(n_treated_work[d] >= n_earlystop){
        break
      }
      
      # Treat one cohort at the current dose level d and generate DLT outcomes
      n_dlt_work[d] = n_dlt_work[d] + rbinom(1, cohort_size, tox_prob_work[d])
      n_treated_work[d] = n_treated_work[d] + cohort_size
      
      #---------------------- SKBD pseudo-posterior update at dose d: ------------------------------------
      #         borrow information from doses with accrued data via kernel-based weights
      
      # indices of doses with observed data (including current dose)
      obs_idx = which(n_treated_work > 0)  
      
      # kernel similarities between current dose value and all dose values on the working grid
      ker_vals = kernel(dose_set_work[d], dose_set_work, symmetric, theta)  # length = n_dose_work
      
      # normalize kernel values over observed doses to form borrowing weights
      ker_sum = sum(ker_vals[obs_idx])
      weight = ker_vals[obs_idx] / ker_sum
      
      # posterior parameters (Beta) after borrowing
      post_alpha_work[d] = pri_alpha + sum(weight * n_dlt_work[obs_idx])
      post_beta_work[d] = pri_beta + sum(weight * (n_treated_work[obs_idx] - n_dlt_work[obs_idx]))

      #------------- determine if the dose insertion conditions are satisfied------------------------------
      res = insert_check(d, post_alpha_work, post_beta_work, key_L, key_U, C1, C2, n_treated_work)
      need_insert = res$need_insert
      insert_code = res$insert_code
      LOC_BELOW_MIN = res$LOC_BELOW_MIN
      LOC_ABOVE_MAX = res$LOC_ABOVE_MAX
      
      #-------------- implement the dose insertion if conditions are satisfied ----------------------------- 
      if (need_insert && sum(n_treated_work) < n_patients && max(n_dlt_work)>0){
        if (insert_code == LOC_BELOW_MIN){
          if (dose_set_work[d] == dose_set_std[1] && sum(n_treated_work) <= (n_patients - 2*cohort_size)){
            tmp <- sort(c(dose_set_work, dose_set_work[d]/2), index.return=T)
            dose_set_work <- tmp$x
            if (n_dose_work==length(tox_prob)){insert_at_cohort[trial]<-sum(n_treated_work)/cohort_size} ## record after how many cohorts the insertion takes place
            
            tox_prob_work <- c(tox_prob_work, ins.p(umi)); tox_prob_work<-tox_prob_work[tmp$ix] 
            n_dose_work <- n_dose_work+1
            
            n_dlt_work <- c(n_dlt_work, 0); n_treated_work <- c(n_treated_work, 0);   
            n_dlt_work <- n_dlt_work[tmp$ix]; n_treated_work <- n_treated_work[tmp$ix]; 
            post_alpha_work <- c(post_alpha_work, pri_alpha); post_beta_work <- c(post_beta_work, pri_beta)
            post_alpha_work <- post_alpha_work[tmp$ix]; post_beta_work <- post_beta_work[tmp$ix]
            d <- 1; is_eliminated<-c(is_eliminated,0); is_eliminated <-is_eliminated[tmp$ix]
            n_insertions_trial <- n_insertions_trial + 1
            ## treat the next cohort of patients with the newly inserted dose
            
            n_dlt_work[d] <- n_dlt_work[d] + rbinom(1,cohort_size, tox_prob_work[d])
            n_treated_work[d] <- n_treated_work[d] + cohort_size
            post_alpha_work[d]<- n_dlt_work[d]+pri_alpha
            post_beta_work[d]<-n_treated_work[d]-n_dlt_work[d]+pri_beta   
            
            bandwidth<-(0.5*abs(dose_set_work[d]-dose_set_work[d+1]))/((1-weight^(1/3))^(1/3))
            fit <- locfit(
              (n_dlt_work+0.05)~lp(dose_set_work,h=bandwidth),weights=n_treated_work+0.5,deg=1,  family="binomial", link="logit",kern="tcub",maxk=5000
            )
            newdata<-seq(dose_set_work[d],dose_set_work[d+1],0.01)###
            dif.new<-abs(pava(predict(fit,newdata))-target)
            newdose<-newdata[min(which(dif.new==min(dif.new)))]
          } else {stop<-1; newdose<-NA}  
        }
        
        
        
        
        
      }
      
      
      
      
      
      
      
      
      
      
      
      
       
      if (need_insert==1 & sum(n_treated_work) < n_patients & max(n_dlt_work)>0){
        if (umi==0){
          if (dose_set_work[d]==dose_set_std[1] & sum(n_treated_work) <= (n_patients-2*cohort_size)){
            tmp <- sort(c(dose_set_work, dose_set_work[d]/2), index.return=T)
            dose_set_work <- tmp$x
            if (n_dose_work==length(tox_prob)){insert_at_cohort[trial]<-sum(n_treated_work)/cohort_size} ## record after how many cohorts the insertion takes place
            
            tox_prob_work <- c(tox_prob_work, ins.p(umi)); tox_prob_work<-tox_prob_work[tmp$ix] 
            n_dose_work <- n_dose_work+1
            
            n_dlt_work <- c(n_dlt_work, 0); n_treated_work <- c(n_treated_work, 0);   
            n_dlt_work <- n_dlt_work[tmp$ix]; n_treated_work <- n_treated_work[tmp$ix]; 
            post_alpha_work <- c(post_alpha_work, pri_alpha); post_beta_work <- c(post_beta_work, pri_beta)
            post_alpha_work <- post_alpha_work[tmp$ix]; post_beta_work <- post_beta_work[tmp$ix]
            d <- 1; is_eliminated<-c(is_eliminated,0); is_eliminated <-is_eliminated[tmp$ix]
            n_insertions_trial <- n_insertions_trial + 1
            ## treat the next cohort of patients with the newly inserted dose
            
            n_dlt_work[d] <- n_dlt_work[d] + rbinom(1,cohort_size, tox_prob_work[d])
            n_treated_work[d] <- n_treated_work[d] + cohort_size
            post_alpha_work[d]<- n_dlt_work[d]+pri_alpha
            post_beta_work[d]<-n_treated_work[d]-n_dlt_work[d]+pri_beta   
            
            bandwidth<-(0.5*abs(dose_set_work[d]-dose_set_work[d+1]))/((1-weight^(1/3))^(1/3))
            fit <- locfit(
              (n_dlt_work+0.05)~lp(dose_set_work,h=bandwidth),weights=n_treated_work+0.5,deg=1,  family="binomial", link="logit",kern="tcub",maxk=5000
            )
            newdata<-seq(dose_set_work[d],dose_set_work[d+1],0.01)###
            dif.new<-abs(pava(predict(fit,newdata))-target)
            newdose<-newdata[min(which(dif.new==min(dif.new)))]
          } else {stop<-1; newdose<-NA}  
        } else if (d!=1 & umi==d-1){
          bandwidth<-(0.5*(dose_set_work[d]-dose_set_work[d-1]))/((1-weight^(1/3))^(1/3))
          fit <- locfit((n_dlt_work+0.05)~lp(dose_set_work, h=bandwidth),weights=n_treated_work+0.5,deg=1,  family="binomial", link="logit",kern="tcub",maxk=5000)
          newdata<-seq(dose_set_work[d-1],dose_set_work[d],0.01)
          dif.new<-abs(pava(predict(fit,newdata))-target)
          newdose<-newdata[min(which(dif.new==min(dif.new)))]
        } else if (umi==d){
          bandwidth<-(0.5*(dose_set_work[d+1]-dose_set_work[d]))/((1-weight^(1/3))^(1/3))
          fit <- locfit((n_dlt_work+0.05)~lp(dose_set_work, h=bandwidth),weights=n_treated_work+0.5,deg=1,  family="binomial", link="logit",kern="tcub",maxk=5000)
          newdata<-seq(dose_set_work[d],dose_set_work[d+1],0.01)
          dif.new<-abs(pava(predict(fit,newdata))-target)
          newdose<-newdata[min(which(dif.new==min(dif.new)))]
        } else {
          if (dose_set_work[d]==1.5*dose[n_dose]){
            stop=1; break
          }else {
            bandwidth<-( 0.5* (1.5*dose[n_dose]-dose_set_work[d]) )/((1-weight^(1/3))^(1/3))
            fit <- locfit((n_dlt_work+0.05)~lp(dose_set_work, h=bandwidth),weights=n_treated_work+0.5,deg=1,  family="binomial", link="logit",kern="tcub",maxk=5000)
            newdata<-seq(dose_set_work[d],1.5*dose[n_dose],0.01)
            dif.new<-abs(pava(predict(fit,newdata))-target)
            newdose<-newdata[min(which(dif.new==min(dif.new)))]
          }
        }
        ## rearrange the dose and toxicity rate vectors to reflect the correct ordering
        if ( (!is.na(newdose)) && (!(newdose %in% dose_set_work)) ){          
          tmp <- sort(c(dose_set_work, newdose), index.return=T)
          dose_set_work <- tmp$x
          if (n_dose_work==length(tox_prob)){insert_at_cohort[trial]<-sum(n_treated_work)/cohort_size} ## record after how many cohorts the insertion takes place
          
          ins.p <- function(umi){ 
            if (d!=1 & umi==d-1){
              if (tox_prob_work[d-1]<target & target<tox_prob_work[d]){
                target
              }else{      
                (tox_prob_work[d-1]+tox_prob_work[d])/2
              }   
            }else if (umi==0){
              if (target<tox_prob_work[d]) {
                target
              }else{
                tox_prob_work[d]/2
              }
            }else if (umi==d){
              if (tox_prob_work[d]<target & target<tox_prob_work[d+1]){
                target
              }else{     
                (tox_prob_work[d]+tox_prob_work[d+1])/2
              }
            }else{
              if (target>tox_prob_work[d]){
                target
              }else{
                tox_prob_work[d]+0.1
              } 
            } 
          }
          
          tox_prob_work <- c(tox_prob_work, ins.p(umi)); tox_prob_work <- tox_prob_work[tmp$ix] 
          n_dose_work <- n_dose_work+1
          
          n_dlt_work <- c(n_dlt_work, 0); n_treated_work <- c(n_treated_work, 0) 
          n_dlt_work <- n_dlt_work[tmp$ix]; n_treated_work <- n_treated_work[tmp$ix]
          
          post_alpha_work <- c(post_alpha_work, pri_alpha); post_beta_work <- c(post_beta_work, pri_beta)
          post_alpha_work <- post_alpha_work[tmp$ix]; post_beta_work <- post_beta_work[tmp$ix]
          is_eliminated<-c(is_eliminated,0); is_eliminated <-is_eliminated[tmp$ix]
          
          d <- which(newdose==dose_set_work)
          if (length(d)>1){
            d <- sample(d,1)
          }
          n_insertions_trial<-n_insertions_trial+1
          ## treat the next cohort of patients with the newly inserted dose
          if(sum(n_treated_work)<n_patients){
            n_dlt_work[d] <- n_dlt_work[d] + rbinom(1,cohort_size, tox_prob_work[d])
            n_treated_work[d] <- n_treated_work[d] + cohort_size
            post_alpha_work[d]<- n_dlt_work[d]+pri_alpha
            post_beta_work[d]<-n_treated_work[d]-n_dlt_work[d]+pri_beta                
          }
        } 
      }
      
      ## BOIN design: determine if the current dose should be eliminated
      if (!is.na(b.elim[n_treated_work[d]])){
        if (n_dlt_work[d]>=b.elim[n_treated_work[d]]){
          is_eliminated[d:n_dose_work]<-1
          if(d==1) {is_earlystop=TRUE; break}
        }
        ## BOIN design: implement the extra safe rule by decreasing the elimination cutoff for the lowest dose
        if (extrasafe){
          if (d==1 && n_dlt_work[1]>3){
            if (1-pbeta(target,n_dlt_work[1]+1,n_treated_work[1]-n_dlt_work[1]+1)>cutoff.eli-offset) {
              is_earlystop=TRUE;  break
            }
          }
        }
      }
      
      ## dose escalation/de-escalation based on BOIN design
      if (n_dlt_work[d]<=b.e[n_treated_work[d]] && d!=n_dose_work) {
        if (is_eliminated[d+1]==0) d<-d+1
      }else if (n_dlt_work[d]>=b.d[n_treated_work[d]] && d!=1) {
        d<-d-1
      }else{d=d}
      
      # stop when total enrolled patients reaches n_patients
      total<-sum(n_treated_work)
      if(total >= n_patients){
        break
      }  
    }#end while
    
    ## dose selection by BOIN design  
    if(is_earlystop) { 
      sel_dose_idx[trial]=99
    }else  {
      sel_dose_idx[trial]=select.mtd(target, n_treated_work, n_dlt_work, cutoff.eli, extrasafe, offset)$MTD 
    }
    ## save the results
    n_insertions[trial]<-n_insertions_trial
    sel_dose[trial]<- ifelse(sel_dose_idx[trial]==0, NA, dose_set_work[sel_dose_idx[trial]])
    temp.simdata <- cbind(rep(trial, length(n_dlt_work)), dose_set_work, n_treated_work, n_dlt_work, rep(sel_dose_idx[trial], length(n_dose_work)))
    simdata <- rbind(simdata, temp.simdata)
  }
  ## output results
  simdata <- data.frame(simdata)
  names(simdata) <- c("Simulation", "Dose", "N", "X", "Selection")
  selpercent<- rep(0, n_dose) # selection percentage at prespecified doses
  ptspercent<- rep(0, n_dose) # percentage of patients at prespecified doses
  for (i in 1:n_dose) {
    selpercent[i]<-sum(sel_dose==dose[i],na.rm=TRUE)/n_trial*100
    ptspercent[i]<-(sum(simdata[which(simdata$Dose==dose[i]),3])/(n_trial*n_patients))*100
  }
  ins.select<- 100-sum(selpercent) # selection percentage at inserted dose
  ins.pts<- 100-sum(ptspercent)    # percentage of patients at inserted dose
  
  ins.dose<-NULL
  for (i in 1:n_trial){
    if (sel_dose[i] %in% dose) {ins.dose<-ins.dose}
    else {ins.dose<-c(ins.dose,sel_dose[i])}
  }
  ins.mean <- mean(ins.dose,na.rm=TRUE)
  ins.sd <- sd(ins.dose, na.rm=TRUE)
  ins.percent <-(1000-sum(n_insertions==0))/n_trial*100 
  cohort.mean <- mean(insert_at_cohort[which(insert_at_cohort!=0)])
  insTimes.median <- median(n_insertions)
  cat ("selection percentage at each prespecified dose level (%):\n")
  cat (formatC(selpercent,digits=2, format="f"),sep="  ","\n")
  cat ("percentage of patients treated at prespecified dose level (%):\n")
  cat (formatC(ptspercent,digits=2, format="f"),sep="  ","\n")
  cat ("mean of the inserted dose (SD):",formatC(ins.mean,digits=2,format="f"),"(",formatC(ins.sd, digits=2, format="f"),") \n")
  cat ("selection percentage of the inserted dose (%):",formatC(ins.select,digits=2, format="f"),"\n")
  cat ("percentage of patients treated at the inserted doses (%):",formatC(ins.pts, digits=2, format="f"),"\n")
  cat ("percentage of trials with dose insertion (%):",formatC(ins.percent,digits=1, format="f"),"\n")
  cat ("The average of cohorts after which trial insertion is likely to take place:",formatC(cohort.mean,digits=1,format="f"),"\n")
  cat ("Median number of insertion times:", formatC(insTimes.median,digits=1,format="f"),"\n")

}