from sklearn.svm import SVC
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.grid_search import GridSearchCV
from glob import glob
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix
outer_folds = 5
inner_folds = 3
path_name = "all"
n_jobs = 5

for index in range(1,1000,2):
#for index in range(2,1002,2):
#for index in range(34,64,1):
#for index in range(64,101,1):
	print index
	pathways = glob('perm_test/permutation_test_{}/'.format(index) + path_name + "_path/*.txt")
	labels = pd.read_csv("tcgaPathways/labels",index_col = 0,sep=" ")
	y = LabelEncoder().fit(labels.x).transform(labels.x)
	n_classes = len(np.unique(y))

	clf = SVC(probability = True,shrinking=False)
	#tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4,1e-5],'C': [1e-3,1e-2,1e-1,1, 10, 100, 1000]},
	#		    {'kernel': ['poly'], 'degree':[1,2,3], 'C': [1e-3,1e-2,1e-1,1, 10, 100, 1000]},
	#		    {'kernel': ['linear'], 'C': [1e-3,1e-2,1e-1,1, 10, 100, 1000]}]

	tuned_parameters = [{'kernel': ['linear'], 'C': [1e-3,1e-2,1e-1,1, 10, 100, 1000]}]

	gs = GridSearchCV(clf,tuned_parameters,n_jobs=n_jobs,cv = inner_folds, refit=True)
	test_predictions = np.zeros((y.shape[0], n_classes*len(pathways)))
	test_scores = {}#np.zeros((len(pathways), outer_folds))
	test_avg_scores = {}

	for i, p in enumerate(pathways):
		print p
		p_name = p.split('.')[0].split('/')[-1]
		X= pd.read_csv(p,index_col = 0,sep=" ")
		test_scores[p_name] = []

		mean_cm = np.zeros((4,4))

		for j, (train, test) in enumerate(StratifiedKFold(y,outer_folds,True)):
			X_train = X.iloc[train,]
			X_test = X.iloc[test,]

			y_train = y[train]
			y_test = y[test]

			gs.fit(X_train,y_train)
			pd.DataFrame(gs.grid_scores_).to_csv('perm_test/permutation_test_{}/'.format(index) + path_name + '_res/' + p_name + '_fold{}_grid_scores.txt'.format(j))
			test_scores[p_name].append(gs.score(X_test,y_test))
			test_predictions[np.ix_(test, range(i * n_classes, (i + 1) * n_classes))] = gs.predict_proba(X_test)
			cm = confusion_matrix(y_test,gs.predict(X_test))
			mean_cm = mean_cm + cm
			np.savetxt('perm_test/permutation_test_{}/'.format(index) + path_name + '_res/' + p_name + '_ranked_genes_fold{}.txt'.format(j), gs.best_estimator_.coef_)
	
		#print mean_cm
		np.savetxt('perm_test/permutation_test_{}/'.format(index) + path_name + '_res/' + p_name + 'test_predictions.txt', test_predictions)
		np.savetxt('perm_test/permutation_test_{}/'.format(index) + path_name + '_res/' + p_name + '_conf_mat.txt', mean_cm)
		test_avg_scores[p_name] = dict(avg_acc=np.mean(test_scores[p_name]),best_params = gs.best_params_)

	pd.DataFrame(test_avg_scores).T.to_csv('perm_test/permutation_test_{}/'.format(index) + path_name + '_res/test_avg_accuracy.txt')
	pd.DataFrame(test_scores).T.to_csv('perm_test/permutation_test_{}/'.format(index) + path_name+ '_res/test_scores.txt')
	pd.DataFrame(test_predictions, index=labels.index).to_csv('perm_test/permutation_test_{}/'.format(index) + path_name+ '_res/level1_features.txt')
