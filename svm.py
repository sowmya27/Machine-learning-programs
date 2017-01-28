__author__ = 'annapurnaannadatha'
import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm, preprocessing
import pandas as pd
from matplotlib import style
style.use("ggplot")
import warnings  # to ignore all the deprecated warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)



data_headers = ['HMM','SS','OGS']

data_df = pd.DataFrame.from_csv("metamorphic.csv")
malware_X = np.array(data_df[data_headers].values)
print(malware_X.shape)
test_mal_category = np.empty([len(malware_X),1])
for i in range(0,len(malware_X)):
    test_mal_category[i]= 1

data1_df = pd.DataFrame.from_csv("benign.csv")
benign_X = np.array(data1_df[data_headers].values)
print(benign_X.shape)
test_ben_category = np.empty([len(benign_X),1])
for i in range(0,len(benign_X)):
    test_ben_category[i]= -1

test_size = 20
train_set = np.concatenate((malware_X[:test_size],benign_X[:test_size]))
test_set = np.concatenate((malware_X[test_size:],benign_X[test_size:]))

train_set = preprocessing.scale(train_set)
test_set = preprocessing.scale(test_set)

train_out = np.concatenate((test_mal_category[:test_size],test_ben_category[:test_size]))
test_out = np.concatenate((test_mal_category[test_size:],test_ben_category[test_size:]))


#Kernel = Polynomial ,P=2,C=0
clf = svm.SVC(kernel="poly",degree =2, C= 0.000001)
clf.fit(train_set,train_out.ravel())
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Kernel = Polynomial learning machine ,P=2,C=0")
print(clf.predict(test_set))

#testing the data and computing the accuracy of the result
correct_count = 0
for x in range(1, len(test_set)+1):
    if clf.predict(test_set[-x]) == test_out[-x]:
        correct_count += 1
print("Misclassification count :", 40 - correct_count )
print("Correct set count :",correct_count)
print("Accuracy :", (correct_count/len(test_set)) * 100.00)
#print(clf.score(train_set,train_out.ravel()))

#Kernel = Polynomial ,P=2,C=3
clf = svm.SVC(kernel="poly",degree =2, C=3)
clf.fit(train_set,train_out.ravel())
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Kernel = Polynomial learning machine,P=2,C=3")
print(clf.predict(test_set))

#testing the data and computing the accuracy of the result
correct_count = 0
for x in range(1, len(test_set)+1):
    if clf.predict(test_set[-x]) == test_out[-x]:
        correct_count += 1
print("Misclassification count :", 40 - correct_count )
print("Correct set count :",correct_count)
print("Accuracy:", (correct_count/len(test_set)) * 100.00)
#print(clf.score(train_set,train_out.ravel()))

#Kernel = radial basis, C = 0
clf = svm.SVC(kernel="rbf",C=0.0001)
clf.fit(train_set,train_out.ravel())
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Kernel = Gaussian radial-basis function ,C=0")
print(clf.predict(test_set))

#testing the data and computing the accuracy of the result
correct_count = 0
for x in range(1, len(test_set)+1):
    if clf.predict(test_set[-x]) == test_out[-x]:
        correct_count += 1
print("Misclassification count :", 40 - correct_count )
print("Correct set count :",correct_count)
print("Accuracy :", (correct_count/len(test_set)) * 100.00)
#print(clf.score(train_set,train_out.ravel()))


#Kernel = radial basis, C = 3
clf = svm.SVC(kernel="rbf",C=3)
clf.fit(train_set,train_out.ravel())
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Kernel = Gaussian radial-basis function,C=3")
print(clf.predict(test_set))

#testing the data and computing the accuracy of the result
correct_count = 0
for x in range(1, len(test_set)+1):
    if clf.predict(test_set[-x]) == test_out[-x]:
        correct_count += 1
print("Misclassification count :", 40 - correct_count )
print("Correct set count :",correct_count)
print("Accuracy using radial basis Kernel,c=3:", (correct_count/len(test_set)) * 100.00)
#print(clf.score(train_set,train_out.ravel()))
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")




'''
w = clf.coef_[0]
a = -w[0] / w[1]
xx = np.linspace(min(test_set[:,0]),max(test_set[:,0]))
yy = a * xx - clf.intercept_[0] / w[1]

h0 =plt.plot(xx,yy,"k-",label="spam/ham")
plt.scatter(test_set[:,0],test_set[:,1],c=test_out)
plt.legend()
plt.show()
'''





