import copy
import pandas as pd
import numpy as np

def RandomSvmModel(self):
    sample = "../Script/Mapping/Mfe/Sample.csv"
    df = pd.read_csv(sample)  # 返回一个DataFrame的对象，这个是pandas的一个数据结构
    cols = list(df.columns.values)
    cols_data = copy.deepcopy(cols)
    # cols_data.remove('sequence')
    cols_data.remove('class')
    x_data = df[list(cols_data)]  # 抽取作为训练数据的各属性值
    x = np.array(x_data)
    print(len(x))
    y_data = df['class']  # 最后一列作为每行对应的标签label
    y = np.array(y_data)
    print(y)
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.25, random_state=33)
    ############################ -*- coding: 得到特征重要性排序 -*-######################
    from sklearn.model_selection import GridSearchCV
    from sklearn.ensemble import RandomForestClassifier
    param_test1 = {'n_estimators': range(100, 1000, 100)}
    gsearch1 = GridSearchCV(estimator=RandomForestClassifier(min_samples_split=100,
                                                             min_samples_leaf=20, max_depth=8, max_features='sqrt',
                                                             random_state=10),
                            param_grid=param_test1, scoring='roc_auc', cv=5)
    gsearch1.fit(X_train, y_train)

    param_test2 = {'max_depth': range(3, 14, 2), 'min_samples_split': range(2, 200, 10)}
    gsearch2 = GridSearchCV(estimator=RandomForestClassifier(n_estimators=350, min_samples_split=10,
                                                             min_samples_leaf=20, max_depth=8, max_features='sqrt',
                                                             random_state=10),
                            param_grid=param_test2, scoring='roc_auc', cv=5)
    gsearch2.fit(X_train, y_train)

    # 再对内部节点再划分所需最小样本数min_samples_split和叶子节点最少样本数min_samples_leaf一起调参
    param_test3 = {'min_samples_split': range(10, 150, 20), 'min_samples_leaf': range(10, 60, 10)}
    gsearch3 = GridSearchCV(estimator=RandomForestClassifier(n_estimators=350, min_samples_split=42,
                                                             min_samples_leaf=20, max_depth=9, max_features='sqrt',
                                                             random_state=10),
                            param_grid=param_test3, scoring='roc_auc', cv=5)
    gsearch3.fit(X_train, y_train)

    param_test4 = {'max_features': range(3, 50, 2)}
    gsearch4 = GridSearchCV(estimator=RandomForestClassifier(n_estimators=350, min_samples_split=10,
                                                             min_samples_leaf=10, max_depth=9, random_state=10),
                            param_grid=param_test4, scoring='roc_auc', cv=5)
    gsearch4.fit(X_train, y_train)
    gsearch4.best_params_
    feat_lables = cols

    forest = RandomForestClassifier(n_estimators=350, min_samples_split=10,
                                    min_samples_leaf=10, max_depth=9, random_state=10, max_features=50,
                                    oob_score=True)
    forest.fit(X_train, y_train)
    importance = forest.feature_importances_
    imp_result = np.argsort(importance)[::-1]
    feature = []
    x_d = []
    y_d = []
    for i in range(X_train.shape[1]):
        temp = imp_result[i]
        # print("%2d. %-*s %f" % (i, 30, cols_data[temp], importance[temp]))
        x_d.append(cols_data[temp])
        y_d.append(importance[temp])
        feature.append(cols_data[temp])
    print(feature[0: 30])
    import matplotlib.pyplot as plt
    plt.bar(x_d[0:30], height=y_d[0:30], tick_label=x_d[0:30])
    plt.setp(plt.gca().get_xticklabels(), rotation=80, horizontalalignment='right')
    plt.show()

    sample = "../Script/Mapping/Mfe/Sample.csv"
    df = pd.read_csv(sample)  # 返回一个DataFrame的对象，这个是pandas的一个数据结构
    cols = list(df.columns.values)
    cols_data = copy.deepcopy(cols)
    # cols_data.remove('sequence')
    cols_data.remove('class')
    x_data = df[list(feature[0: 49])]  # 抽取作为训练数据的各属性值
    x = np.array(x_data)
    y_data = df['class']  # 最后一列作为每行对应的标签label
    y = np.array(y_data)
    import matplotlib.pyplot as plt
    plt.bar(range(len(importance[0:30])), importance[0:30], tick_label=cols_data[0:30])
    plt.setp(plt.gca().get_xticklabels(), rotation=80, horizontalalignment='right')
    plt.show()

    from sklearn import svm
    from sklearn.model_selection import GridSearchCV
    import matplotlib.pyplot as plt

    from sklearn.metrics import roc_curve
    from sklearn.metrics import roc_auc_score
    from sklearn.metrics import classification_report

    param_test1 = {
        'C': [x for x in np.arange(0.5, 2, 0.1)],
        'gamma': [x for x in np.arange(0.01, 0.1, 0.01)],
        'tol': [x for x in np.arange(0.01, 0.1, 0.01)],
    }

    scores = ['precision', 'recall']

    print("# Tuning hyper-parameters for roc_auc")
    print()
    grid = GridSearchCV(estimator=svm.SVC(probability=True),
                        param_grid=param_test1, scoring='accuracy', refit=True, cv=5)
    # accuracy,roc_auc
    grid.fit(X_train, y_train)
    print(grid.best_params_)
    print(grid.best_score_)
    print("Grid scores on development set:")
    print()
    means = grid.cv_results_['mean_test_score']
    stds = grid.cv_results_['std_test_score']
    for mean, std, params in zip(means, stds, grid.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r" % (mean, std * 2, params))

    print()
    print("Detailed classification report:")
    print()
    print("The model is trained on the full development set.")
    print("The scores are computed on the full evaluation set.")
    print()
    y_true, y_pred = y_test, grid.predict(X_test)
    # 打印在测试集上的预测结果与真实值的分数
    print(classification_report(y_true, y_pred))
    print("================================")
    print(classification_report(y_train, grid.predict(X_train)))

    print(grid.best_estimator_)
    y_pred_pro = grid.predict_proba(X_test)
    y_scores = pd.DataFrame(y_pred_pro, columns=grid.classes_.tolist())[1].values
    auc_value = roc_auc_score(y_true, y_scores)
    # clf = svm.SVC(kernel='rbf', class_weight='balanced')
    # clf.fit(X_train, y_train)
    # 绘制ROC曲线
    fpr, tpr, thresholds = roc_curve(y_true, y_scores, pos_label=1.0)
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='blue', label='AUC = %0.4f' % auc_value)
    plt.plot([0, 1], [0, 1], color='red', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")
    plt.show()