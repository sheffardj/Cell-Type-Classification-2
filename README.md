# Cell-Type-Classification-2 (FINAL PROJECT)
For Univeristy of Victoria Course: Stat 454 | Machine Learning 

This was an iteration on /Cell-Type-Classification-1 (MIDTERM PROJECT) with a harder to classify dataset, and we we're restricted to 10 five-fold Cross Validations. We started with a Random Forest and XGBoost, but the Random Forest did not compete, so we created an Ensemble Method with the XGBoost and an Elastic Net Regression, with Elastic Net predictions replacing low-confidence XGBoost results. 

We performed a class-label pair-wise feature reduction. We plotted the F1-scores of XGBoost/Elastic-Net (XGB-EN) versus just the Elastic Net Regression (EN) (which was the best performer of the /Cell-Type-Classification-1 project). Finally, we created a 6 page report in academic style. 
