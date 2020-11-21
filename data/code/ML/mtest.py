#coding=utf-8
# import tensorflow as tf
import tensorflow.compat.v1 as tf

tf.disable_v2_behavior()
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


def MLPdemo():
    # 数据格式处理
    # sample = "../Script/Mapping/Mfe/Sample.csv"
    sample = "../Script/Mapping/Mfe/SamTr.csv"
    potus = list(csv.reader(open(sample)))
    dx = []; dy = []
    potus = potus[1:]
    # shuffle(potus)
    for i in range(0, len(potus)):
        dx.append([int(x) for x in potus[i][0:len(potus[i]) - 1]])
        dy.append([int(potus[i][len(potus[i]) - 1])])
    # train_dx = dx[0:864]; test_dx = dx[1152:]
    train_dx = dx[0:864]; test_dx = dx[864:]
    train_dy = dy[0:864]; test_dy = dy[864:]


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
        # print(str(w1.eval(session=sess)))
        # print(str(w2.eval(session=sess)))
        # print(b1.eval(session=sess))
        # print(b2.eval(session=sess))
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
        print("auc  :", roc_test,"-", roc_train)
        y_scores = prob_sig
        fpr, tpr, thresholds = roc_curve(test_dy, prob, pos_label=1.0)
        plt.figure(figsize=(6.4, 6.4))
        plt.plot(fpr, tpr, color='blue', label='AUC = %0.4f' % roc_test)
        plt.plot([0, 1], [0, 1], color='red', linestyle='--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver operating characteristic of MLP')
        plt.legend(loc="lower right")
        plt.show()

if __name__ == "__main__":
    MLPdemo()
