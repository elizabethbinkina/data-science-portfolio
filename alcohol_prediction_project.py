# -*- coding: utf-8 -*-
"""alcohol_prediction_project.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1KyFI8QrnD7nYUeAh9jFhpJW-Z_j8kSAL

# What Variables Influence Alchohol Consumption

Elizabeth Binkina

## 1. Introduction
"""

#imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import bernoulli
from scipy.stats import binom
from scipy.stats import norm
from scipy import stats
from scipy.stats import t

import statsmodels.api as sm
import statsmodels.formula.api as smf

import seaborn as sns; sns.set()
from scipy.stats import f

import math

"""# Motivation:

Student alcohol consumption is influenced by many factors, and influences others. In this dataset, we are trying to explore the explanatory variables that contribute the most to student drinking, as well as what variables could be affected. A dataset of this type can be helpful for people to understand what exactly is the most powerful influence of student drinking, as well as what kind of consequences result from it. For example, “ About one in four college students report experiencing academic difficulties from drinking, such as missing class or getting behind in schoolwork” (USDHHS 3). In this dataset, we are trying to figure out if this is true, as well as if all the other factors that can contribute to consuming more alcohol.

“College Drinking.” National Institute on Alcohol Abuse and Alcoholism, U.S. Department of Health and Human Services, https://www.niaaa.nih.gov/publications/brochures-and-fact-sheets/college-drinking.

# Research Questions
1. Is there a difference in number of absences when consuming more alcohol on the weekends for students who want a higher education?

We plan to answer this question by making a scatter plot with the line of best fit to explore the relationship between weekend drinking and absences for students who want a higher education vs those who dont. This will be interesting to find out if wanting a higher education affects both attendance and the amount of alcohol consumption over the weekend.

1. Is there an association between extra curricular activities and alcohol consumption during the weekdays?

We plan to answer the question by doing a hypothesis test for difference in means to answer whether there is an association between being involved in extra curricular activity and workday alcohol consumption. We would want to know if there is an association because many students participate in after school activities and it would be interesting to know if that affects the amount they drink during the week.


2. Is there a linear relationship between students' grade averages and weekend alcohol consumption, workday alcohol consumption, free time after school, and health status?

We plan to answer this question by making a linear regression model and equation to answer which variable have the strongest linear relationship with grade average. In other words, we will find out what variables affect overall grade averages in students the most. We would be interested to know what variable affects grade averages among students the most.


4. What explanatory variables should we include in the model that predicts weekend alcohol consumption to build a parsimonious model?

We plan to answer this question by building a logistic regression model using the Backwards Elimination Algorithm to see what factors are most prominent in increased alcohol consumption in students.

5. Is there a linear relationship between the log-odds of the success level of weekend alcohol consumption and final student grades , going out with friends, quality of family relationships, and amount of time spent studying in the math student dataset?

We plan to answer this question to first conducting a hypothesis test to see if the slopes are non-zero or not in the logistic regression model. Then, find if there is a linear relationship with the log-odds of the success level of weekend alcohol consumption.

# Dataset :
This data set is a random survey of math students in secondary schools with social, gender, and study information. It contains information about there family, their parents jobs, freetime, social life, their health, freetime, their grades, and how much alcohol they consume on the weekdays and weekends.
"""

df = pd.read_csv('student-mat.csv')
df

"""-------------------------------------------------------------

## 2. Descriptive Analytics

# Research Question : Is there a difference in number of absences when consuming more alcohol on the weekends for students who want a higher education?

These responce variables are chosen because we are going to be able to figure out if students who consume more alcohol during the weekend skip class more often, as well as if wanting a higher education affects the amount of alcohol consumption during the weekend.
"""

df['Walc'].hist()
plt.title('Weekend Alcohol Consumption')
plt.ylabel('Frequency')
plt.xlabel('Weekend Alcohol Consumption')
plt.show()

"""Alcohol consumption on weekends tends to be skewed right. This implies that more people drink less on the weekends."""

df['absences'].hist()
plt.title('Number of Absences')
plt.ylabel('Frequency')
plt.xlabel('Absences')
plt.show()

"""Number of absences is skewed right. This implies that most students in this dataset go to class more often than not."""

df[['Walc', 'absences']].corr()

"""Drinking on the weekends and number of absences do not have a strong linear relationship."""

df_alc_counts = df['Walc'].value_counts()
df_alc_counts

"""Measure of Center:"""

df_alc_counts.median()

df['Walc'].median()

"""The median level of drinking for students is 2 out of 5, which represents about 80 students.

Measure of Spread:
"""

alc_iqr = df_alc_counts.quantile(0.75)-df_alc_counts.quantile(0.25)
alc_iqr

"""Because the distrobution is skewed right, the median and IQR are more appropraite measure of center and spread to use.

The IQR is 34.
"""

plot1 = sns.lmplot(x= 'Walc', y= 'absences', hue = 'higher', data = df)
plt.title('Weekend Alcohol Consumption vs Number of Absences')
plt.xlabel("Weekend Alcohol Consumption")
plt.show()

"""Students who consume more alcohol on the weekends and want a higher education, tend to have more absences than students who consume more alcohol on the weekends and do not want a higher education. The intercept of the best fit line of students who want a higher education is about 4, meaning that students who do not drink on the weekends and want a higher education, have on average 4 absences. The slope is about 0.8, meaning if student's level of alcohol consumed on the weekends increases by 1, their number of absences will increase by 0.8 on average. The intercept of the best fit line of students who do not want a higher education is about 13, meaning that students who do not drink on the weekends and do not want a higher education, have on average of 13 absences. The slope is about -2.4, meaning that if student's level of alcohol consumption increased by 1,their number of absences will decrease by 2.4 on average.
In this model, there exists 3 large outliars. Since the dataset has just under 400 values, it will not skew the line of best fit greatly. Therfore, there is a difference in number of absences in students based on whether they want a higher education and how much they drink on the weekends.

-----------------------------------------------

## 3. Inference
"""

df_yes = df[df.activities == "yes"]
df_no = df[df.activities == "no"]

"""# Research Question: Is there sufficient evidence to suggest that there is an association between extra curricular activity involvement and workday alcohol consumption?

Extra curricular involvement is being tested because it is often promoted that the more involved in extra curricular activities you are, the less likely you are to get involved with bad activities like drinking. Workday consumption was chosen because that is the main topic being tested in this study.

H0: u1 - u2 = 0

HA: u1 - u2 ≠ 0

u1 = average workday alcohol consumption for students involved in extra curricular activities

u2 = average workday alcohol consumption for students not involved in extra curricular activites

Hypothesis Test Conditions:

The first condition is met because the sample sizes of both dataframes are greater than 30 at 201 and 194.

The second condition is met because the samples were both randomly sample.

The third condition is met because both samples sizes are less than 10% of the population.

The fourth condition is met because the samples are independent of each other.

# Hypothesis Test:
"""

mean_yes = df_yes.Dalc.mean()
mean_no = df_no.Dalc.mean()
std_yes = df_yes.Dalc.mean()
std_no = df_yes.Dalc.mean()

se = (((std_yes**2)/201) + ((std_no**2)/194))**.5
se

t1 = (mean_yes - mean_no - 0) / se
t1

p_value = t.cdf(t1, 393)
p_value

"""Conclusion:

Because the p-value of .205 is greater than our alpha of .05, we fail to reject the null hypothesis. Therefore, we don't have sufficient evidence that there is an association between extra curricular activity involvement and workday alcohol consumption.

----------------------------

## 4. Linear Regression
"""

df["Ga"] = (df.G1 + df.G2 + df.G3) / 3

"""# Research Question: Is there a linear relationship between grade average and weekend alcohol consumption, workday alcohol consumption, free time after school, and health status?

I chose the two alcohol consumption variables because I mainly want to test the affect that alcohol consumption has on grade average. I did, however, include two other variables, free time and health status, to see if there was another variable that had a strong impact on grade average.

List of Variables Used:
* x = Ga
* y = walc
* dalc
* freetime
* health
"""

results = smf.ols("Ga ~ Walc + Dalc + freetime + health", data=df).fit()
results.summary()

"""# Linear Regression Equation:

Grade Average-hat = 11.6896 - 0.1891(Weekend Alcohol Consumption) - 0.1222(Workday Alcohol Consumption) + .0935(Free Time) - .1966(Health)
"""

plt.scatter(results.fittedvalues, results.resid)
plt.ylabel("Residual")
plt.xlabel("Fitted Value")

plt.hist(results.resid)
plt.xlabel("Residuals")
plt.show()

sns.pairplot(df)

"""# Linear Regression Conditions:

The linearity condition is met because the residuals are clumped around y=0.

The constant variance of residuals condition is met because residuals are distributed fairly evenly accross.

The residuals are normal condition isn't met because although the graph looks symmetrical, it is not normal due to a drop in residuals around 0.

The independence condition is met because they were randomly sampled and the sample size is less than 10% of the population.

The non colinearity condition is also met when looking at all the graphs.

# Fit of the Model:

R-Squared Method: The r-squared of this model is .014 which shows that the model does not fit very well. Only 1.4% of the variance is explained by the model.

Linearity Method: When looking at the residuals graph, the residuals are centered pretty well at y=0. This shows that the model has a very good fit.

None of the slopes in the model have sufficient evidence to suggest that they are non-zero. This is because when looking at the summary output, the p-values testing if each of the variables are non-zero are all greater than our alpha of .05. Therefore, there is not sufficient evidence for each variable to prove that each variable is non-zero.

# Conclusion:

Research Question: "Is there a linear relationship between grade average and weekend alcohol consumption, workday alcohol consumption, free time after school, and health status?"

Based on the two methods testing the fit of the model, we can not come to a conclusion on the answer of this question. When utilizing the two methods, the r-squared method showed there was non fit, but the linearity method showed there was a fit. Therefore, we can not come to a conclusion on whether there is a linear relationship in this model.

----------------------------

## 5. Logistic Regression

# Research Question: What explanatory variables should we include in the model that predicts weekend alcohol consumption to build a parsimonious model?

I chose the response variable to be weekend alcohol consumption, because I'm interested to see what factors are most prominent in increased alcohol consumption on students. I declined to use weekday alcohol consumption, because I feel that number would be a poor indicator of how many students actually use alcohol, as I believe most students are most likley to use alcohol on the weekends, rather than the weekdays. Since weekend alcohol consumption was a numerical variable with 1 being the lowest consumed and 5 being the highest, I chose to classify people who put 3 and above as people who used alcohol, and below 3 as people who didn't. I chose 3 as my threshold because I wanted to classify people who used alcohol both moderately and heavily as people who used alcohol in my model. I chose not to include students who recorded a 2 or lower in the group consuming alcohol, because I feel students who don't drink at all or have a few drinks on the weekends are more responsible, and prioritize activities and studies over partying.

I chose G3(Final Grade, numeric), goout(Going out with friends, numeric), age(Student's age, numeric), famrel(Quality of family relationships, numeric), and study time(weekly study time, numeric) as my explanatory variables. I wanted to see if these variables have any effect on how much a student drinks, specifically how school work, performance, home lifes, and life outside of school affect alcohol consumption. I chose G3 and studytime variables to see if lower grades and less study time played into increased alcohol usage. I chose age to see if there was a strong relationship between older students and more alcohol use. I chose famrel to see if poor family relationships resulted in higher usage of alcohol. Finally, I chose the goout variable to see if students who go out more also choose to drink more alcohol.
"""

df.columns

"""I first needed to create a 0/1 categorical variable from my numeric response variable, WALC, so I created it using the threshold of 3."""

df['y'] = 1*(df['Walc']>=3)

"""I then separated my dataset into training and test datasets to train and then test our logistic regression model."""

from sklearn.model_selection import train_test_split
df_train, df_test = train_test_split(df, test_size = 0.2)

len(df_test)/len(df)

"""I began implementing a Backwards Elimination Algorithm to find the most parsimonious model, and see which combination of the explanatory variables I chose will result in the most parsimonious model."""

current_mod = smf.logit('y~G3+goout+age+famrel+studytime', data=df_train).fit()
print('Iteration 1: AIC of Current Model', current_mod.aic)

# Model without G3
test_mod = smf.logit('y~goout+age+famrel+studytime', data=df_train).fit()
print('AIC of Test Model that Deletes G3 from Current Model', test_mod.aic)

# Model without goout
test_mod = smf.logit('y~G3+age+famrel+studytime', data=df_train).fit()
print('AIC of Test Model that Deletes goout from Current Model', test_mod.aic)

# Model without age
test_mod = smf.logit('y~G3+goout+famrel+studytime', data=df_train).fit()
print('AIC of Test Model that Deletes age from Current Model', test_mod.aic)

# Model without famrel
test_mod = smf.logit('y~G3+goout+age+studytime', data=df_train).fit()
print('AIC of Test Model that Deletes famrel from Current Model', test_mod.aic)

# Model without studytime
test_mod = smf.logit('y~G3+goout+age+famrel', data=df_train).fit()
print('AIC of Test Model that Deletes studytime from Current Model', test_mod.aic)

"""The model without age had the lowest AIC, so we continue our backwards elimination process."""

current_mod = smf.logit('y~G3+goout+famrel+studytime', data=df_train).fit()
print('Iteration 2: AIC of Current Model', current_mod.aic)

# Model without G3
test_mod = smf.logit('y~goout+famrel+studytime', data=df_train).fit()
print('AIC of Test Model that Deletes G3 from Current Model', test_mod.aic)

# Model without goout
test_mod = smf.logit('y~G3+famrel+studytime', data=df_train).fit()
print('AIC of Test Model that Deletes goout from Current Model', test_mod.aic)

# Model without famrel
test_mod = smf.logit('y~G3+goout+studytime', data=df_train).fit()
print('AIC of Test Model that Deletes famrel from Current Model', test_mod.aic)

# Model without studytime
test_mod = smf.logit('y~G3+goout+famrel', data=df_train).fit()
print('AIC of Test Model that Deletes studytime from Current Model', test_mod.aic)

"""We stop the backwards elimination process, because no other model had a lower AIC then the current model. We then display our final logistic regression model."""

final_mod = smf.logit('y~G3+goout+famrel+studytime', data=df_train).fit()
final_mod.summary()

"""Final Logistic Regression Equation:

$log(\frac{\hat{p}}{1-\hat{p}}) = -1.4806 + 0.0521G3 + 0.8308goout - 0.3063famrel - 0.5069studytime$

We will next calculate the ROC Curve and AUC using the logistic regression curve on our test dataset. First, we must retrieve the predicitive probabilities of the test dataset with the trained model.
"""

phat_test = final_mod.predict(exog=df_test)
phat_test.head()

df_test['phat_test']=phat_test
df_test.head()

from sklearn.metrics import confusion_matrix, roc_curve, roc_auc_score
fprs, tprs, thresholds = roc_curve(y_true=df_test['y'],
                            y_score=df_test['phat_test'])
auc = roc_auc_score(y_true=df_test['y'],
                    y_score=df_test['phat_test'])
print(auc)

def plot_roc(fpr, tpr, auc, lw=2):
    plt.plot(fpr, tpr, color='darkorange', lw=lw,
             label='ROC curve (area = '+str(round(auc,3))+')')
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend(loc="lower right")
    plt.show()

plot_roc(fprs, tprs, auc)

"""The graph shows there is a predictive probability threshold that comes close to the top left of the graph, and this model will provide accurate classifications of the test dataset. The AUC is strong at 0.824, meaning there is a predictive probability threshold that comes close to the ideal false positive, true positive ratio of (0,1). We will next use the ROC to find a good predictive probability threshold. Based on the research question, I'd like a predictive probability threshold that correctly identifies most true positives, without focusing too much on false positives. Using the fact that 0.5 is the most common threshold, I will go a little bit lower and put it at 0.4, because we want to maximize true positives.  """

df_test['yhat'] = 1*(df_test['phat_test']>0.4)
df_test.head()

tn, fp, fn, tp = confusion_matrix(y_true=df_test['y'],
                                  y_pred=df_test['yhat']).ravel()
(tn, fp, fn, tp)

"""Based on the threshold of 0.4, we were given tpr of 0.82 and an fpr of 0.32.

# Research Question: Is there a linear relationship between the log-odds of the success level of weekend alcohol consumption and final student grades , going out with friends, quality of family relationships, and amount of time spent studying in the math student dataset?

In order to answer the above question, we must conduct a hypothesis test on whether the slopes are non-zero or not in the logistic regression model. First, we must check conditions for conducting the hypothesis test.

Condition 1: Independence of sample
The sample is randomly collected, and the n of 395<10% of all math students

Condition 2: Linearity Condition
All curves have some semblance of an S shape, so this condition is met
"""

sns.lmplot(x="G3", y='y',data=df, logistic=True)
plt.ylabel('1=approve, 0=disapprove')
plt.show()

sns.lmplot(x="goout", y='y',data=df, logistic=True)
plt.ylabel('1=approve, 0=disapprove')
plt.show()

sns.lmplot(x="studytime", y='y',data=df, logistic=True)
plt.ylabel('1=approve, 0=disapprove')
plt.show()

sns.lmplot(x="famrel", y='y',data=df, logistic=True)
plt.ylabel('1=approve, 0=disapprove')
plt.show()

"""Condition 3: No Multi-Collinearity
This condition is not met because the explanatory variables all have a linear relationship. We will conduct hypothesis tests anyway
"""

sns.pairplot(df[['G3', 'goout', 'famrel', 'studytime']])

"""H0 = G3_slope = goout_slope = famrel_slope = studytime_slope = 0
HA = G3_slope, goout_slope, famrel_slope, studytime_slope != 0
"""

final_mod.summary()

"""Using an alpha value equal to 0.05, we have reson to reject the null hypothesis for the slopes of goout, famrel, and studytime. We have sufficient evidence to suggest the alternative, that the slopes of goout, famrel, and studytime are not 0. However, we fail to reject the null hypothesis on the G3 slope, and we do not have sufficient evidence to suggest the slope of G3 is not 0.

To answer the first of our two research questions, we found that the slopes of G3, goout, famrel, and studytime should be the explanatory variables to build a parsimonious model that predicts weekend alcohol consumption. For the second question, we found the variables of goout, famrel, and studytime have a linear relationship with the log-odds of the success level of weekend alcohol consumption, while the variable G3 does not.

## 6. Conclusion

# Summary:
 We have learned a lot from our research into the math student dataset, and we were able to find answers to all of our research questions. While attempting to identify whether there was a difference in absences for people who drank on the weekends while also seeking higher education, we found a significant difference. Students who consume alcohol but want a higher education have more absences, while those who don't want higher education have fewer absences. This suggests students are more stressed while drinking and seeking higher education, resulting in more absences. In the search to find an association between extra-curricular activities and weekday alcohol consumption, our hypothesis test did not conclude there was one. This finding was interesting, as we were expecting participation in extra-curriculars to curb weekday drinking, but found that not to be the case. In an effort to identify a linear relationship between grade average and other variables, we weren't able to correctly identify one. Grade average doesn't have a strong linear connection with either weekday or weekend drinking, amount of freetime, or health status, suggesting we may need to search for other variables to find a relationship. Finally, we asked if there was a relationship between the log-odds success of weekend alcohol consumption and other explanatory variables, and found the variables of going out with friends, amount of time spent studying, and family relationships all had a significant impact. Going back to our motivation, we found that wanting a higher education and excess consumption of alcohol can lead to more absences, suggesting the stress that comes from higher schooling can lead to more drinking. We also found extra curriculars don't have as much of an impact as one might think to reduce drinking, and that good family relationships, less studying, and more likelihood to go out with friends will increase alcohol intake.

# Future Work:
To continue our discussion on student alcohol consumption, we could ask questions such as "Is there a linear relationship with alcohol consumption and parental factors? We could ask if parental factors have any affect on if students drink alcohol during the weekend or weekdays. We might also want to analyze different cultures/geographic regions and how they affect student alcohol consumption.
Shortcomings: Some shortcomings include the fact that the dataset was a little small, and there was some ambiguity as to where the dataset was collected. I'd like to see a more wide-ranging dataset that covers numerous different areas, so we can get better analysis and get expected results in our dataset.
"""