/******************************************************************/
/*  STUDENT PROJECT: DESCRIPTIVE AND INFERENTIAL ANALYSIS         */
/*  PROJECT: Diabetes Risk Prediction                              */
/*  Group: 5                                                       */
/*  FILES: Diabetes_S2.csv (data set for odd-numbered teams)       */
/******************************************************************/

/*-----------------------*/
/*  PART 1: DESCRIPTIVE  ANALYSIS   */
/*-----------------------*/

/* Step 1: Import the dataset */
PROC IMPORT DATAFILE="/home/u64136089/Mycode/Final proj/diabetes_S2.csv"  
            OUT=diabetes
            DBMS=CSV 
            REPLACE;
    GETNAMES=YES;  /* Use the first row as variable names */
RUN;

/* Step 2: Descriptive Statistics for Continuous Variables */
/* We find the dataset includes the following continuous variables:
   BMI MentHlth PhysHlth */
PROC MEANS DATA=diabetes N MEAN STD MIN MAX;
    VAR BMI MentHlth PhysHlth;
    TITLE "Descriptive Statistics for Continuous Variables";
RUN;

/*Plot histogram for continuous variables.*/
proc sgplot data=diabetes;
	histogram BMI;
	title "BMI distribution";
run;

proc sgplot DATA=diabetes;
	histogram MentHlth;
	title "Mental Health Days Distribution";
run;

proc sgplot DATA=diabetes;
	histogram PhysHlth;
	title "Physical Health Days Distribution";
run;

/* Step 3: Frequency Analysis for Categorical Variables */
PROC FREQ DATA=diabetes;
    TABLES Sex Fruits Veggies Smoker Stroke HighBP HighChol CholCheck  DiffWalk
           NoDocbcCost PhysActivity HvyAlcoholConsump AnyHealthcare 
           Diabetes_binary HeartDiseaseorAttack GenHlth Age Education Income / 
           NOCUM NOROW NOCOL;
    TITLE "Frequency Tables for Key Categorical Variables";
RUN;

/* Plot bar Chart for Categorical Variables */
%MACRO bar_chart(var);
   PROC SGPLOT DATA=diabetes;
      VBAR &var;
      TITLE "Bar Chart for &var";
   RUN;
%MEND;

%bar_chart(Diabetes_binary);
%bar_chart(HighBP);
%bar_chart(HighChol);
%bar_chart(CholCheck);
%bar_chart(Smoker);
%bar_chart(Stroke);
%bar_chart(HeartDiseaseorAttack);
%bar_chart(PhysActivity);
%bar_chart(Fruits);
%bar_chart(Veggies);
%bar_chart(HvyAlcoholConsump);
%bar_chart(AnyHealthcare);
%bar_chart(NoDocbcCost);
%bar_chart(DiffWalk);
%bar_chart(Sex);
%bar_chart(Age);
%bar_chart(Education);
%bar_chart(Income);
%bar_chart(GenHlth);


/*-------------------------*/
/* PART 2: INFERENTIAL ANALYSIS    */
/*-------------------------*/

/* Step 4 : Analysis for Continuous Variables */
/* Since BMI is approximately normally distributed, we'll use t-tests */
PROC TTEST DATA=diabetes;
    CLASS Diabetes_binary;
    VAR BMI;
    TITLE "T-test for BMI by Diabetes Status";
RUN;

/* For MentHlth and PhysHlth which are likely non-normal, use Wilcoxon test */
PROC NPAR1WAY DATA=diabetes WILCOXON;
    CLASS Diabetes_binary;
    VAR MentHlth PhysHlth;
    TITLE "Wilcoxon Rank Sum Tests for Mental and Physical Health";
RUN;

/* Step 5: Analysis for Binary Categorical Variables */
/* Use chi-square tests for binary variables */
%MACRO chi_test(var);
    PROC FREQ DATA=diabetes;
        TABLES Diabetes_binary*&var / CHISQ;
        TITLE "Association between Diabetes and &var";
    RUN;
%MEND;

%chi_test(HighBP);
%chi_test(HighChol);
%chi_test(CholCheck);
%chi_test(Smoker);
%chi_test(Stroke);
%chi_test(HeartDiseaseorAttack);
%chi_test(PhysActivity);
%chi_test(Fruits);
%chi_test(Veggies);
%chi_test(HvyAlcoholConsump);
%chi_test(AnyHealthcare);
%chi_test(NoDocbcCost);
%chi_test(DiffWalk);
%chi_test(Sex);

/* Step 6: Analysis for Ordinal/Multi-level Categorical Variables */
/* For ordinal variables (GenHlth, Age, Education, Income), use Cochran-Armitage trend test */
PROC FREQ DATA=diabetes;
    TABLES Diabetes_binary*GenHlth / TREND;
    TABLES Diabetes_binary*Age / TREND;
    TABLES Diabetes_binary*Education / TREND;
    TABLES Diabetes_binary*Income / TREND;
    TITLE "Trend Tests for Ordinal Variables";
RUN;

/* For variables where expected cell counts <5, use Fisher's exact test */
PROC FREQ DATA=diabetes;
    TABLES Diabetes_binary*HvyAlcoholConsump / FISHER;
    TABLES Diabetes_binary*Stroke / FISHER;
    TITLE "Fisher's Exact Tests for Variables with Small Cell Counts";
RUN;

/* Step 7: Fit the Best Model using Logistic Regression with Variable Selection */
proc logistic data=diabetes descending;
   class HighBP HighChol CholCheck Smoker Stroke HeartDiseaseorAttack 
         PhysActivity Fruits Veggies HvyAlcoholConsump AnyHealthcare 
         NoDocbcCost DiffWalk Sex GenHlth Age Education Income / param=ref;
   model Diabetes_binary = HighBP HighChol CholCheck Smoker Stroke 
         HeartDiseaseorAttack PhysActivity Fruits Veggies HvyAlcoholConsump 
         AnyHealthcare NoDocbcCost DiffWalk Sex GenHlth Age Education Income 
         BMI MentHlth PhysHlth
         / selection=stepwise slentry=0.05 slstay=0.05 details;
   title "Stepwise Logistic Regression to Predict Diabetes Status";
run;

/*-----------------------*/
/*  PART 3                */
/*-----------------------*/

/* Step 8: Fit the GLM */
proc genmod data=diabetes descending;
   class GenHlth(ref='5') 
         HighBP(ref='0') 
         HighChol(ref='0') 
         DiffWalk(ref='0') 
         Age(ref='13') / param=ref;
   model Diabetes_binary = GenHlth HighBP HighChol DiffWalk Age BMI
         / dist=binomial link=logit;
   title "GLM Model: Reference Level for Age = 13 (oldest group)";
run;

/*Step 9: Multinomial Logistic Regression (for GenHlth and Age)*/
options locale=en_US;
proc logistic data=diabetes descending;
   class GenHlth (ref='5') Age (ref='13') / param=ref;
   model Diabetes_binary = GenHlth Age;
   oddsratio GenHlth;
   oddsratio Age;
   title "Logistic Regression with Ordinal Predictors: GenHlth and Age";
run;

/*Step 10: Two-Way Interactions Terms and Effect Modification*/
proc logistic data=diabetes descending;
   class GenHlth HighBP HighChol DiffWalk Age / param=ref;
   model Diabetes_binary = 
      GenHlth HighBP HighChol DiffWalk Age BMI
      GenHlth*HighBP GenHlth*HighChol GenHlth*DiffWalk GenHlth*Age GenHlth*BMI
      HighBP*HighChol HighBP*DiffWalk HighBP*Age HighBP*BMI
      HighChol*DiffWalk HighChol*Age HighChol*BMI
      DiffWalk*Age DiffWalk*BMI
      Age*BMI
      / selection=stepwise slentry=0.05 slstay=0.05;
   title "Main Effects + 2-Way Interactions with Stepwise Selection";
run;

/*Step 11: Principal Component Analysis*/
/* Step 11.1: Prepare numeric-only dataset for PCA */
data pca_data;
   set diabetes;
   if DiffWalk = 1 then WalkDiff_num = 1;
   else if DiffWalk = 0 then WalkDiff_num = 0;
run;

/* Step 11.2: Principal Component Analysis */
proc princomp data=pca_data out=pca_scores std;
   var BMI Age WalkDiff_num;
   title "PCA on BMI, Age, DiffWalk (numeric)";
run;

/* Step 11.3: Logistic regression using 1st two principal components */
proc logistic data=pca_scores descending;
   model Diabetes_binary = Prin1 Prin2;
   title "Logistic Regression using Principal Components";
run;

/* Step 12: Model validation - ROC, AUC, Classification stats */
proc logistic data=diabetes descending plots(only)=roc;
   class GenHlth(ref='5') HighBP(ref='0') HighChol(ref='0') DiffWalk(ref='0') Age(ref='13') / param=ref;
   model Diabetes_binary = GenHlth HighBP HighChol DiffWalk Age BMI;
   score data=diabetes out=pred_data;
   title "Logistic Model Validation: ROC Curve and Classification Stats";
run;

/* Step 13: Classification table and predicted probabilities */
proc logistic data=diabetes descending;
   class GenHlth(ref='5') HighBP(ref='0') HighChol(ref='0') DiffWalk(ref='0') Age(ref='13') / param=ref;
   model Diabetes_binary = GenHlth HighBP HighChol DiffWalk Age BMI;
   output out=val_results p=pred_prob;
run;

proc freq data=val_results;
   tables Diabetes_binary*pred_prob / nopercent norow nocol;
run;

/*------------------------------------------------------*/
/* PART 4: Additional Analyses to Enhance Results       */
/*------------------------------------------------------*/

/* Step 14: ROC Curve and AUC to Evaluate Final Model */
proc logistic data=diabetes descending plots(only)=roc(id=prob);
   class HighBP HighChol DiffWalk GenHlth Age / param=ref;
   model Diabetes_binary = HighBP HighChol DiffWalk GenHlth Age BMI;
   roc;
   title "ROC Curve and AUC for Final Logistic Model Performance";
run;

/* Step 15: Confusion Matrix - Classification Table */
proc logistic data=diabetes descending;
   class HighBP HighChol DiffWalk GenHlth Age / param=ref;
   model Diabetes_binary = HighBP HighChol DiffWalk GenHlth Age BMI;
   output out=predicted p=Pred_Prob predprobs=individual;
   title "Model Scoring: Creating Predicted Probabilities (with _INTO_)";
run;

proc freq data=predicted;
   tables Diabetes_binary*_INTO_ / norow nocol nopercent;
   title "Confusion Matrix: Actual vs Predicted Diabetes Status";
run;

/* Step 16: Influence Diagnostics - Identify Outliers */
proc logistic data=diabetes descending;
   class HighBP HighChol DiffWalk GenHlth Age / param=ref;
   model Diabetes_binary = HighBP HighChol DiffWalk GenHlth Age BMI;
   output out=diag_results pred=pred p=prob xbeta=logit resdev=devres h=leverage;
   title "Influence Diagnostics for Final Logistic Model";
run;
proc sgplot data=diag_results;
   scatter x=leverage y=devres / markerattrs=(symbol=circlefilled);
   refline 0 / axis=y;
   title "Diagnostic Plot: Deviance Residuals vs Leverage";
run;
proc means data=diag_results max;
   var leverage devres;
run;

/* Step 17: Check Multicollinearity - Variance Inflation Factor (VIF) */
proc reg data=diabetes;
   model Diabetes_binary = HighBP HighChol DiffWalk GenHlth Age BMI / vif;
   title "VIF Check: Multicollinearity Assessment Among Predictors";
run;

/* Step 18: Hosmer-Lemeshow Goodness-of-Fit Test */
proc logistic data=diabetes descending;
   class HighBP HighChol DiffWalk GenHlth Age / param=ref;
   model Diabetes_binary = HighBP HighChol DiffWalk GenHlth Age BMI / lackfit;
   title "Hosmer-Lemeshow Goodness-of-Fit Test for Final Model";
run;

/* Step 19: Sensitivity and Specificity Analysis at Cutoff = 0.5 */
data predicted;
    set predicted;
    if Pred_Prob >= 0.5 then PredictedClass = 1;
    else PredictedClass = 0;
run;

proc freq data=predicted;
    tables Diabetes_binary*PredictedClass / nopercent norow nocol;
    title "Sensitivity and Specificity at 0.5 Cutoff";
run;

