from sklearn.svm import SVC
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.grid_search import GridSearchCV
from glob import glob
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import LabelEncoder

outer_folds = 5
inner_folds = 3
path_name = "all78"
n_jobs = 2

labels = pd.read_csv("labels_gbm.txt",index_col = 0,sep=" ")
X = pd.read_csv(path_name + '_res/level1_features.txt',index_col = 0)
y = LabelEncoder().fit(labels.x).transform(labels.x)
n_classes = len(np.unique(y))

clf = SVC(probability = True,shrinking=False)
#tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4,1e-5],'C': [1e-3,1e-2,1e-1,1, 10, 100, 1000]},
#		    {'kernel': ['poly'], 'degree':[1,2,3], 'C': [1e-3,1e-2,1e-1,1, 10, 100, 1000]},
#		    {'kernel': ['linear'], 'C': [1e-3,1e-2,1e-1,1, 10, 100, 1000]}]

tuned_parameters = [{'kernel': ['linear'], 'C': [1e-3,1e-2,1e-1,1, 10, 100, 1000]}]

gs = GridSearchCV(clf,tuned_parameters,n_jobs=n_jobs,cv = inner_folds, refit=True)
test_scores = []
test_avg_scores = {}


for i, (train, test) in enumerate(StratifiedKFold(y,outer_folds,True)):
	X_train = X.iloc[train,]
	X_test = X.iloc[test,]

	y_train = y[train]
	y_test = y[test]

	gs.fit(X_train,y_train)
	print i
	pd.DataFrame(gs.grid_scores_).to_csv(path_name + '_res/level2_fold{}_grid_scores.txt'.format(i))
	test_scores.append(gs.score(X_test, y_test))
	np.savetxt(path_name + '_res/ranked_pathways_fold{}.txt'.format(i), gs.best_estimator_.coef_)

np.savetxt(path_name + '_res/level2_test_accuracies.txt', test_scores)
