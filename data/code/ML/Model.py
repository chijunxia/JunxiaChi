#coding=utf-8
import csv
import copy
import logging
import matplotlib.pyplot as plt
from ML.feature import feature, cols_data, details, redundancy, rmdup, cols_data_rmdup
import numpy as np
import pandas as pd
from sklearn.metrics import classification_report
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import recall_score, accuracy_score, f1_score, precision_score

logging.getLogger().setLevel(logging.INFO)


import csv
import logging
from random import shuffle
import os
logging.getLogger().setLevel(logging.INFO)
# os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
# os.environ["CUDA_VISIBLE_DEVICES"] = '4' #use GPU with ID=0
# config = tf.ConfigProto()
# config.gpu_options.per_process_gpu_memory_fraction = 0.5 # maximun alloc gpu50% of MEM
# config.gpu_options.allow_growth = True #allocate dynamically


class Model:
    def __init__(self):
        self.count = 0
        self.Dicts = {}

    @staticmethod
    def DataInit(evalute=1):
        from sklearn.preprocessing import StandardScaler
        from sklearn.model_selection import train_test_split
        if evalute == 0:
            df = pd.read_csv("../Script/Mapping/Mfe/SamTrain.csv")
            cols = list(df.columns.values)
            cols_data = copy.deepcopy(cols)
            cols_data.remove('class')
            x_data = df[list(cols_data)]
            x = np.array(x_data)
            x = StandardScaler().fit_transform(x)
            y_data = df['class']
            y = np.array(y_data)
            data = np.insert(x, x[0].size, values=y, axis=1)
            # np.random.shuffle(data)
            y = data[:, data[0].size - 1]
            x = np.delete(data, -1, axis=1)
            X_train, X_test, Y_train, Y_test = train_test_split(x, y, test_size=0.25, random_state=1)
        elif evalute==1:
            # dftrain = pd.read_csv("../Script/Mapping/Mfe/SamTr.csv")
            dftrain = pd.read_csv("../Script/Mapping/Mfe/AllSampleTrain_0314.csv")
            cols = list(dftrain.columns.values)
            cols_data = copy.deepcopy(cols)
            cols_data.remove('class')
            cols_data.remove('sequence')
            # for i in rmdup:
            #     cols_data.remove(i)
            x_data = dftrain[list(cols_data)]
            x = np.array(x_data)
            # X_train = x
            X_train = StandardScaler().fit_transform(x)
            y_data = dftrain['class']
            Y_train = np.array(y_data)

            # dftest = pd.read_csv("../Script/Mapping/Mfe/SamTe.csv")
            dftest = pd.read_csv("../Script/Mapping/Mfe/AllSampleTest_0314.csv")
            cols = list(dftest.columns.values)
            cols_data = copy.deepcopy(cols)
            cols_data.remove('class')
            cols_data.remove('sequence')
            # for i in rmdup:
            #     cols_data.remove(i)
            x_data = dftest[list(cols_data)]
            x = np.array(x_data)
            # X_test = x
            X_test = StandardScaler().fit_transform(x)
            y_data = dftest['class']
            Y_test = np.array(y_data)
        elif evalute== 2:
            df = pd.read_csv("../Script/Mapping/Mfe/SampleR.csv")
            cols = list(df.columns.values)
            cols_data = copy.deepcopy(cols)
            cols_data.remove('class')
            x_data = df[list(cols_data)]
            x = np.array(x_data)
            x = StandardScaler().fit_transform(x)
            y_data = df['class']
            y = np.array(y_data)
            data = np.insert(x, x[0].size, values=y, axis=1)
            # np.random.shuffle(data)
            y = data[:, data[0].size - 1]
            x = np.delete(data, -1, axis=1)
            X_train, X_test, Y_train, Y_test = train_test_split(x, y, test_size=0.25, random_state=1)
            data_train = np.insert(X_train, X_train[0].size, values=Y_train, axis=1)
            data_test = np.insert(X_test, X_test[0].size, values=Y_test, axis=1)
            # with open("../Script/Mapping/Mfe/SamTr.csv", "w") as f:
            #     f_csv = csv.writer(f)
            #     f_csv.writerow(cols)
            #     f_csv.writerows(data_train)
            # with open("../Script/Mapping/Mfe/SamTe.csv", "w") as ft:
            #     f_csv = csv.writer(ft)
            #     f_csv.writerow(cols)
            #     f_csv.writerows(data_test)
            # return
        else:
            return "输入不合法"

        return X_train, X_test, Y_train, Y_test

    @staticmethod
    def FigResult(mod, X_train, Y_train, X_test, Y_test, m_name):
        print(m_name)
        Y_test_predict = mod.predict(X_test)
        Y_train_predict = mod.predict(X_train)
        acc_scores = cross_val_score(mod, X_train, Y_train, cv=5, scoring='accuracy')
        rec_scores = cross_val_score(mod, X_train, Y_train, cv=5, scoring='recall')
        fon_scores = cross_val_score(mod, X_train, Y_train, cv=5, scoring='f1')
        pre_scores = cross_val_score(mod, X_train, Y_train, cv=5, scoring='precision')
        # return [list(acc_scores), list(rec_scores), list(pre_scores), list(fon_scores)]
        # print("十折交叉平均值(accuracy) ：", acc_scores.mean())
        # print("十折交叉平均值(recall)   ：", rec_scores.mean())
        # print("十折交叉平均值(f1)       ：", fon_scores.mean())
        # print("十折交叉平均值(precision)：", pre_scores.mean())
        print("accuracy ：", accuracy_score(Y_test, Y_test_predict),"-", accuracy_score(Y_train, Y_train_predict))
        print("recall   ：", recall_score(Y_test, Y_test_predict),"-", recall_score(Y_train, Y_train_predict))
        print("f1       ：", f1_score(Y_test, Y_test_predict),"-", f1_score(Y_train, Y_train_predict))
        print("precision：", precision_score(Y_test, Y_test_predict),"-", precision_score(Y_train, Y_train_predict))

        Y_test_predict_pro = mod.predict_proba(X_test)
        y_scores = pd.DataFrame(Y_test_predict_pro, columns=mod.classes_.tolist())[1].values
        auc_value = roc_auc_score(Y_test, y_scores)
        Y_train_predict_pro = mod.predict_proba(X_train)
        y_train_scores = pd.DataFrame(Y_train_predict_pro, columns=mod.classes_.tolist())[1].values
        auc_value_train = roc_auc_score(Y_train, y_train_scores)
        print("auc  :", auc_value,"-", auc_value_train)
        # print(classification_report(Y_test, Y_test_predict))
        fpr, tpr, thresholds = roc_curve(Y_test, y_scores, pos_label=1.0)
        plt.figure(figsize=(6.4, 6.4))
        plt.plot(fpr, tpr, color='blue', label='AUC = %0.4f' % auc_value)
        plt.plot([0, 1], [0, 1], color='red', linestyle='--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver operating characteristic of '+ m_name)
        plt.legend(loc="lower right")
        plt.show()

    @staticmethod
    def FigAllImg(m1, m2, m3, m4, m5, m6, m7, f, t, s, X_train, Y_train, X_test, Y_test):
        Y_tpp_1 = m1.predict_proba(X_test)
        y_s_1 = pd.DataFrame(Y_tpp_1, columns=m1.classes_.tolist())[1].values
        auc_1 = roc_auc_score(Y_test, y_s_1)
        print(auc_1)
        fpr1, tpr1, thresholds1 = roc_curve(Y_test, y_s_1, pos_label=1.0)

        Y_tpp_2 = m2.predict_proba(X_test)
        y_s_2 = pd.DataFrame(Y_tpp_2, columns=m2.classes_.tolist())[1].values
        auc_2 = roc_auc_score(Y_test, y_s_2)
        print(auc_2)
        fpr2, tpr2, thresholds2 = roc_curve(Y_test, y_s_2, pos_label=1.0)

        Y_tpp_3 = m3.predict_proba(X_test)
        y_s_3 = pd.DataFrame(Y_tpp_3, columns=m3.classes_.tolist())[1].values
        auc_3 = roc_auc_score(Y_test, y_s_3)
        print(auc_3)
        fpr3, tpr3, thresholds3 = roc_curve(Y_test, y_s_3, pos_label=1.0)

        Y_tpp_4 = m4.predict_proba(X_test)
        y_s_4 = pd.DataFrame(Y_tpp_4, columns=m4.classes_.tolist())[1].values
        auc_4 = roc_auc_score(Y_test, y_s_4)
        print(auc_4)
        fpr4, tpr4, thresholds4 = roc_curve(Y_test, y_s_4, pos_label=1.0)

        Y_tpp_5 = m5.predict_proba(X_test)
        y_s_5 = pd.DataFrame(Y_tpp_5, columns=m5.classes_.tolist())[1].values
        auc_5 = roc_auc_score(Y_test, y_s_5)
        print(auc_5)
        fpr5, tpr5, thresholds5 = roc_curve(Y_test, y_s_5, pos_label=1.0)

        Y_tpp_6 = m6.predict_proba(X_test)
        y_s_6 = pd.DataFrame(Y_tpp_6, columns=m6.classes_.tolist())[1].values
        auc_6 = roc_auc_score(Y_test, y_s_6)
        print(auc_6)
        fpr6, tpr6, thresholds6 = roc_curve(Y_test, y_s_6, pos_label=1.0)

        Y_tpp_7 = m7.predict_proba(X_test)
        y_s_7 = pd.DataFrame(Y_tpp_7, columns=m7.classes_.tolist())[1].values
        auc_7 = roc_auc_score(Y_test, y_s_7)
        print(auc_7)
        fpr7, tpr7, thresholds7 = roc_curve(Y_test, y_s_7, pos_label=1.0)


        plt.figure(figsize=(8, 8))
        # plt.plot(fpr1, tpr1, 'r-', label='KNN AUC = %0.4f' % auc_1)
        # plt.plot(fpr2, tpr2, 'k-', label='BAY AUC = %0.4f' % auc_2)
        # plt.plot(fpr3, tpr3, 'y-', label='DTR AUC = %0.4f' % auc_3)
        # plt.plot(fpr4, tpr4, 'g-', label='LRG AUC = %0.4f' % auc_4)
        # plt.plot(fpr5, tpr5, 'c-', label='RFT AUC = %0.4f' % auc_5)
        # plt.plot(fpr6, tpr6, 'b-', label='SVM AUC = %0.4f' % auc_6)
        # plt.plot(fpr7, tpr7, 'm-', label='XGB AUC = %0.4f' % auc_7)
        # plt.plot(f, t, color="orange", label='MLP AUC = %0.4f' % s)

        plt.plot(fpr1, tpr1, 'r-', label='KNN')
        plt.plot(fpr2, tpr2, 'k-', label='BAY')
        plt.plot(fpr3, tpr3, 'y-', label='DTR')
        plt.plot(fpr4, tpr4, 'g-', label='LRG')
        plt.plot(fpr5, tpr5, 'c-', label='RFT')
        plt.plot(fpr6, tpr6, 'b-', label='SVM')
        plt.plot(fpr7, tpr7, 'm-', label='XGB')
        plt.plot(f, t, color="orange", label='MLP')

        plt.plot([0, 1], [0, 1], color='r', linestyle='--')
        plt.plot([0, 0], [0, 1], color='r', linestyle='--')
        plt.plot([1, 0], [1, 1], color='r', linestyle='--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver operating characteristic of All Model')
        plt.legend(loc="lower right")
        plt.show()

    def KnnDemo(self, X_train, X_test, Y_train, Y_test):
        from sklearn.neighbors import KNeighborsClassifier
        mod = KNeighborsClassifier(n_neighbors=57)
        mod.fit(X_train, Y_train)
        Model.FigResult(mod, X_train, Y_train, X_test, Y_test, "KNN")
        return mod

    def Bysdemo(self, X_train, X_test, Y_train, Y_test):
        # GaussianNB 参数只有一个：先验概率priors; 参数的意义主要参考https://www.cnblogs.com/pinard/p/6074222.html
        from sklearn.naive_bayes import GaussianNB, MultinomialNB, BernoulliNB
        mod = GaussianNB()
        mod.fit(X_train, Y_train)
        Model.FigResult(mod, X_train, Y_train, X_test, Y_test, "Bays")
        return mod

        # option = [GaussianNB, BernoulliNB]
        # evalution = ['accuracy', 'recall', 'precision', 'f1']
        #
        # scores = [[], [], [], []]
        # count = -1
        #
        # for e in evalution:
        #     count += 1
        #     for o in option:
        #         mod = o()
        #         mod.fit(X_train, Y_train)
        #         scores[count].append(np.mean(cross_val_score(mod, X_train, Y_train, cv=5, scoring=e)))
        #
        # print(scores)
        # plt.xlabel('Evaluation Criteria')
        # plt.ylabel('Value')
        # plt.ylim(0, 1)
        # num_list = [0.7673668360624883, 0.7739463601532567, 0.7693773688785859, 0.7714057714669738]
        # x = list(range(len(evalution)))
        # total_width, n = 0.8, 3
        # width = total_width / n
        # plt.bar(x, num_list, width=width, label='Multinomial', fc='b')
        # for i in range(len(x)):
        #     x[i] += width
        # num_list2 = [x[0] for x in scores]
        # plt.bar(x, num_list2, width=width, label='Gaussian', tick_label=evalution, fc='g')
        # for i in range(len(x)):
        #     x[i] += width
        # num_list3 = [x[1] for x in scores]
        # plt.bar(x, num_list3, width=width, label='Bernoulli', fc='y')
        # plt.legend()
        # plt.show()


    def Dtree(self, X_train, X_test, Y_train, Y_test):
        from sklearn import tree
        mod = tree.DecisionTreeClassifier(criterion='entropy', max_depth=4)
        mod = mod.fit(X_train, Y_train)
        Model.FigResult(mod, X_train, Y_train, X_test, Y_test, "Dtree")
        return mod

        # scores_en = []
        # scores_gn = []
        # for i in range(1, 50):
        #     mod = tree.DecisionTreeClassifier(criterion='entropy', max_depth=i)
        #     mod.fit(X_train, Y_train)
        #     scores_en.append(np.mean(cross_val_score(mod, X_test, Y_test, cv=3, scoring='accuracy')))
        # print(scores_en)
        # plt.plot(range(len(scores_en)), scores_en, label='entropy', marker='^')
        # for i in range(1, 50):
        #     mod = tree.DecisionTreeClassifier(criterion='gini', max_depth=i)
        #     mod.fit(X_train, Y_train)
        #     scores_gn.append(np.mean(cross_val_score(mod, X_test, Y_test, cv=3, scoring='accuracy')))
        # print(scores_gn)
        # plt.plot(range(len(scores_gn)), scores_gn, label='gini', marker='.', color='r')
        # plt.xlabel("Number of experiments")
        # plt.legend(loc=4)
        # plt.show()
        # return

    def LogisticRegression(self, X_train, X_test, Y_train, Y_test):
        from sklearn.linear_model import LogisticRegression as LR
        mod = LR(penalty="l1", solver="liblinear")
        mod.fit(X_train, Y_train)
        Model.FigResult(mod, X_train, Y_train, X_test, Y_test, "LR")
        return mod

        # lrl1 = LR(penalty="l1", solver="liblinear", C=0.5, max_iter=1000)
        # lrl2 = LR(penalty="l2", solver="liblinear", C=0.5, max_iter=1000)
        #
        # lrl1 = lrl1.fit(X_train, Y_train)
        # lrl2 = lrl2.fit(X_train, Y_train)
        #
        # (lrl1.coef_ != 0).sum(axis=1)  # l1正则化，参数被稀疏化，可以防止过拟合
        # (lrl2.coef_ != 0).sum(axis=1)  # l2正则化，非零稀疏比较多
        # # 调参数
        # l1test = []
        # l2test = []
        # l1 = []
        # l2 = []
        #
        # for i in np.linspace(0.05, 2, 39):
        #     lrl1 = LR(penalty="l1", solver="liblinear", C=i, max_iter=1000)
        #     lrl2 = LR(penalty="l2", solver="liblinear", C=i, max_iter=1000)
        #
        #     lrl1 = lrl1.fit(X_train, Y_train)
        #     l1.append(accuracy_score(lrl1.predict(X_train), Y_train))
        #     l1test.append(accuracy_score(lrl1.predict(X_test), Y_test))
        #
        #     lrl2 = lrl2.fit(X_train, Y_train)
        #     l2.append(accuracy_score(lrl2.predict(X_train), Y_train))
        #     l2test.append(accuracy_score(lrl2.predict(X_test), Y_test))
        #
        # print(l1test)
        # print(l2test)
        # graph = [l1, l2, l1test, l2test]
        # color = ["darkblue", "darkgreen", "deepskyblue", "turquoise"]
        # label = ["L1 train", "L2 train", "L1 test", "L2 test"]
        #
        # for i in range(len(graph)):
        #     plt.plot(np.linspace(0.05, 2, 39), graph[i], color[i], label=label[i])
        # plt.xlabel("Regularization strength")
        # plt.legend()  # loc = 4, 表示右下角
        # plt.show()

    def RFClassdemo(self, X_train, X_test, Y_train, Y_test):
        from sklearn.ensemble import RandomForestClassifier
        mod = RandomForestClassifier(n_estimators=91, max_depth=6, max_features=5)
        mod.fit(X_train, Y_train)
        Model.FigResult(mod, X_train, Y_train, X_test, Y_test, "Random Forest")
        return mod

        from sklearn.model_selection import GridSearchCV
        param_test = {
            'n_estimators': range(1, 102, 10),
            'max_depth': range(3, 14, 2),
            'max_features': range(3, 11, 2),
        }
        mod = GridSearchCV(estimator=RandomForestClassifier(), param_grid=param_test, scoring='accuracy', refit=True, cv=3)
        mod.fit(X_train, Y_train)
        print(mod.best_params_)
        print(mod.best_score_)
        return

    def SVMDemo(self, X_train, X_test, Y_train, Y_test):
        from sklearn import svm
        from sklearn.model_selection import GridSearchCV, learning_curve
        grid = svm.SVC(probability=True, C=0.25, gamma=0.01)
        grid.fit(X_train, Y_train)
        Model.FigResult(grid, X_train, Y_train, X_test, Y_test, "SVM")
        return grid

        param_test = {
            'C': [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            'gamma': [0.001, 0.005, 0.01, 0.02, 0.03, 0.04],
        }
        scores = ['precision', 'recall']
        grid = GridSearchCV(estimator=svm.SVC(probability=True), param_grid=param_test, scoring='accuracy', refit=True, cv=3)
        # grid = svm.SVC(probability=True, C=0.6, gamma=0.01)
        grid.fit(X_train, Y_train)
        print(grid.best_params_)
        print(grid.best_score_)
        return

    # 特征重要性排序——随机森林
    def RFdemo(self):
        X_train, X_test, Y_train, Y_test = Model.DataInit(evalute=1)

        from sklearn.ensemble import RandomForestRegressor
        rf = RandomForestRegressor(n_estimators=81, max_depth=11, max_features=5)
        rf.fit(X_train, Y_train)
        print(sorted(zip(map(lambda x:round(x, 4), rf.feature_importances_), cols_data), reverse=True))
        from matplotlib import pyplot
        pyplot.bar(range(len(rf.feature_importances_)), rf.feature_importances_)
        pyplot.show()

        pyplot.figure(figsize=(12, 6))
        top30 = sorted(zip(map(lambda x:round(x, 4), rf.feature_importances_), cols_data_rmdup), reverse=True)
        top30_title = []
        top30_conte = []
        top30 = top30[0:60]
        for i in top30:
            top30_title.append(i[1])
            top30_conte.append(i[0])
        pyplot.bar(top30_title, top30_conte, tick_label=top30_title)
        plt.setp(plt.gca().get_xticklabels(), rotation=75, horizontalalignment='right')
        pyplot.show()


    # 特征重要性排序——XGBoost
    def XGBoostdemo(self, X_train, X_test, Y_train, Y_test):
        import os
        os.environ["CUDA_VISIBLE_DEVICES"] = '0'
        from xgboost import XGBClassifier
        from matplotlib import pyplot

        # model = XGBClassifier(gpu_id=0, n_estimators=70, max_depth=5, subsample=1)
        model = XGBClassifier(gpu_id=0, n_estimators=90, max_depth=2, subsample=0.7, gamma=0.001)
        model.fit(X_train, Y_train, eval_metric='auc')
        # Model.FigResult(model, X_train, Y_train, X_test, Y_test)

        top_sort = sorted(zip(map(lambda x: round(x, 4), model.feature_importances_), cols_data_rmdup), reverse=True)

        # 画图
        top30 = sorted(zip(map(lambda x: round(x, 4), model.feature_importances_), cols_data_rmdup), reverse=True)[0:30]
        top30_title = []
        top30_conte = []
        map_feature_details = {}
        for i in range(len(feature)):
            map_feature_details[feature[i]] = details[i]
        count = 0
        for i in top30:
            count += 1
            top30_title.append(map_feature_details[i[1]])
            top30_conte.append(i[0])
            print(count, ":", i[1])
        plt.figure(figsize=(8, 6))
        pyplot.bar(top30_title, top30_conte, tick_label=top30_title)
        plt.setp(plt.gca().get_xticklabels(), rotation=90, horizontalalignment='right')
        plt.ylabel("Feature importance ranking")
        pyplot.show()
        return


        X_train_sort = []
        X_test_sort = []
        for k,v in top_sort:
            X_train_sort.append(X_train[:, cols_data.index(v)])
            X_test_sort.append(X_test[:, cols_data.index(v)])

        X_train_sort = np.array(X_train_sort)
        X_train_sort = X_train_sort.transpose()
        X_test_sort = np.array(X_test_sort)
        X_test_sort = X_test_sort.transpose()

        # X_train_i = np.delete(X_train_sort, range(202, 203), axis=1)
        # X_test_i = np.delete(X_test_sort, range(202, 203), axis=1)
        #
        # print(X_train_sort.shape)
        # print(X_test_sort.shape)
        # print(X_train_i.shape)
        # print(X_test_i.shape)
        # Model.XGBoost(X_train_sort, Y_train, X_test_sort, Y_test)
        # return

        parm = []
        t = 0
        max_count = 0
        max_index = 0
        length = 180
        for i in range(1, length):
            X_train_i = np.delete(X_train_sort, range(i, length), axis=1)
            X_test_i = np.delete(X_test_sort, range(i, length), axis=1)
            t = Model.XGBoost(X_train_i, Y_train, X_test_i, Y_test)
            print(t)
            parm.append(t)
            if t > max_count:
                max_count = t
                max_index = i
        print(max_index)
        plt.plot(range(1, len(parm)+1), parm)
        plt.xlabel("Number of features")
        plt.ylabel('accuracy')
        plt.show()


    def XGBoostparm(self, X_train, X_test, Y_train, Y_test):
        from xgboost import XGBClassifier
        import os
        os.environ["CUDA_VISIBLE_DEVICES"] = '0'
        from xgboost import XGBClassifier
        from sklearn.model_selection import GridSearchCV
        mod = XGBClassifier(gpu_id=0, n_estimators=90, max_depth=2, subsample=0.7, gamma=0.001)
        mod.fit(X_train, Y_train)
        Model.FigResult(mod, X_train, Y_train, X_test, Y_test, "XGBoost")
        return mod

        param_test = {
            # 'n_estimators': range(10, 101, 10),
            # 'max_depth': range(3, 12, 1),
            # 'subsample': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        }

        mod = GridSearchCV(estimator=XGBClassifier(gpu_id=0), param_grid=param_test, scoring='accuracy', refit=True, cv=3)
        mod.fit(X_train, Y_train)
        print(mod.best_params_)
        print(mod.best_score_)
        return

        m = XGBClassifier(gpu_id=0)
        m.fit(X_train, Y_train, eval_metric='auc')
        acc_scores = cross_val_score(m, X_train, Y_train, cv=5, scoring='accuracy')
        Model.FigResult(m, X_train, Y_train, X_test, Y_test)
        return acc_scores.mean()

    @staticmethod
    def XGBoost(X_train, Y_train, X_test, Y_test):
        from xgboost import XGBClassifier
        mod = XGBClassifier(gpu_id=0, n_estimators=90, max_depth=2, subsample=0.7, gamma=0.001)
        # mod = XGBClassifier(gpu_id=0, n_estimators=90, max_depth=8, subsample=0.9)
        mod.fit(X_train, Y_train)
        # acc_scores = cross_val_score(mod, X_train, Y_train, cv=3, scoring='accuracy')
        Model.FigResult(mod, X_train, Y_train, X_test, Y_test, "XGBoost")
        # return acc_scores.mean()

    @staticmethod
    def MLPdemo():
        import tensorflow.compat.v1 as tf
        tf.disable_v2_behavior()
        # 数据格式处理
        sample = "../Script/Mapping/Mfe/SamTr.csv"
        potus = list(csv.reader(open(sample)))
        dx = []
        dy = []
        potus = potus[1:]
        # shuffle(potus)
        for i in range(0, len(potus)):
            dx.append([int(x) for x in potus[i][0:len(potus[i]) - 1]])
            dy.append([int(potus[i][len(potus[i]) - 1])])
        train_dx = dx[0:864]
        test_dx = dx[864:]
        train_dy = dy[0:864]
        test_dy = dy[864:]

        # 定义输入和输出
        x = tf.placeholder(tf.float32, shape=(None, 203), name="x-input")
        y_ = tf.placeholder(tf.float32, shape=(None, 1), name="y-input")
        # 定义神经网络的参数
        w1 = tf.Variable(tf.random_normal([203, 10], mean=0, stddev=1, seed=1))
        w2 = tf.Variable(tf.random_normal([10, 1], mean=0, stddev=1, seed=1))
        b1 = tf.Variable(tf.random_normal([10], mean=0, stddev=1, seed=1))
        b2 = tf.Variable(tf.random_normal([1], mean=0, stddev=1, seed=1))

        y1 = tf.matmul(x, w1) + b1
        y11 = tf.nn.relu(y1)
        y2 = tf.matmul(y11, w2) + b2
        y = tf.sigmoid(y2)
        # tf.clip_by_value(t, clip_value_min, clip_value_max,name=None)
        # cross_entropy = -tf.reduce_mean(y_ * tf.log(tf.clip_by_value(y, 1e-10, 1.0)))
        # loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(labels=y_, logits=y))
        loss = -tf.reduce_mean(
            y_ * tf.log(tf.clip_by_value(y, 1e-10, 1.0)) + (1 - y_) * tf.log(tf.clip_by_value(1 - y, 1e-10, 1.0)))
        train_step = tf.train.AdamOptimizer(0.001).minimize(loss)
        X = train_dx
        Y = train_dy
        # 创建会话运行TensorFlow程序
        with tf.Session() as sess:
            init = tf.initialize_all_variables()
            saver = tf.train.Saver()
            sess.run(init)
            steps = 1500
            for i in range(steps):
                # 通过选取样本训练神经网络并更新参数
                sess.run(train_step, feed_dict={x: X, y_: Y})
                # 每迭代1000次输出一次日志信息
                if i % 100 == 0:
                    # 计算所有数据的交叉熵
                    total_cross_entropy, prob = sess.run([loss, y], feed_dict={x: test_dx, y_: test_dy})
                    # 输出交叉熵之和
                    print("After %d training step(s),cross entropy on all data is %g" % (i, total_cross_entropy))
            prob_train = sess.run(y, feed_dict={x: train_dx, y_: train_dy})
            from sklearn.metrics import roc_auc_score, accuracy_score, f1_score, recall_score, precision_score
            roc_test = roc_auc_score(test_dy, prob)
            roc_train = roc_auc_score(train_dy, prob_train)
            prob_sig = []
            for i in prob:
                prob_sig.append(1 if float(i) > 0.5 else 0)
            print(accuracy_score(test_dy, prob_sig))

            # save_path = saver.save(sess, '../ML/model.ckpt')
            # print("Model saved in file: %s" % save_path)
            result = []
            result.append([roc_test, str(w1.eval(session=sess)), str(w2.eval(session=sess)), str(b1.eval(session=sess)),
                           str(b2.eval(session=sess))])

            import matplotlib.pyplot as plt
            from sklearn.metrics import roc_curve
            import pandas as pd
            fpr, tpr, thresholds = roc_curve(test_dy, prob, pos_label=1.0)

            print("auc  :", roc_test,"-", roc_train)
            y_scores = prob_sig
            plt.figure(figsize=(6.4, 6.4))
            plt.plot(fpr, tpr, color='blue', label='AUC = %0.4f' % roc_test)
            plt.plot([0, 1], [0, 1], color='red', linestyle='--')
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('Receiver operating characteristic of MLP')
            plt.legend(loc="lower right")
            plt.show()

            return fpr, tpr, roc_test


if __name__ == '__main__':
    modle = Model()
    # 特征提取
    # modle.preProcess('../Script/Mapping/NavSample_seq.fasta')

    # 数据分类
    X_train, X_test, Y_train, Y_test = Model.DataInit(evalute=1)

    # # 模型选择
    # m1 = modle.KnnDemo(X_train, X_test, Y_train, Y_test)
    # m2 = modle.Bysdemo(X_train, X_test, Y_train, Y_test)
    # m3 = modle.Dtree(X_train, X_test, Y_train, Y_test)
    # m4 = modle.LogisticRegression(X_train, X_test, Y_train, Y_test)
    # m5 = modle.RFClassdemo(X_train, X_test, Y_train, Y_test)
    # m6 = modle.SVMDemo(X_train, X_test, Y_train, Y_test)
    # m7 = modle.XGBoostparm(X_train, X_test, Y_train, Y_test)
    # f, t, s = Model.MLPdemo()
    # Model.FigAllImg(m1=m1, m2=m2, m3=m3, m4=m4, m5=m5, m6=m6, m7=m7, f=f, t=t, s=s, X_train=X_train, X_test=X_test, Y_train=Y_train, Y_test=Y_test)

    # 特征选择
    # modle.RFdemo()
    modle.XGBoostdemo(X_train, X_test, Y_train, Y_test)
    # modle.XGBoost(X_train, Y_train, X_test, Y_test)

