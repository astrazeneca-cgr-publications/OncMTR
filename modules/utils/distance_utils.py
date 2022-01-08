import numpy as np
from scipy.spatial.distance import cdist

	
def cross_entropy(a, b, epsilon=1e-12):
	a = np.array(a)
	b = np.array(b)
	
	#a = np.clip(a, epsilon, 1. - epsilon)
	N = a.shape[0]
	ce = -np.sum(b*np.log(a+1e-9))/N
	
	return ce



def get_cdist_metrics(signal_A, signal_B, ndigits=4):

	signal_A = np.array(signal_A).reshape(1, -1)
	signal_B = np.array(signal_B).reshape(1, -1)
	
	
	distance_metrics = ['euclidean','cityblock','sqeuclidean','cosine','correlation','jaccard','chebyshev','canberra','braycurtis']

	dist_dict = {}

	for dist in distance_metrics:
		dist_dict[dist] = cdist(signal_A, signal_B, dist)[0][0]

	for k, _ in dist_dict.items():
		dist_dict[k] = round(dist_dict[k], ndigits)

	return dist_dict
