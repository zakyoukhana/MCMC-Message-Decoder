import csv
import random
import numpy as np
import os
import sys
import string
import matplotlib.pyplot as plt


class Cipher(object):
	'''
	Substitution Cipher Object

	Designed for use in the Metropolis Hastings Decoder 
	to solve substitution ciphers. Class can be used to encrypt strings,
	randomly generate and transition ciphers, and translate an encrypted
	string using its key
	'''
	def __init__(self, alphabet):
		'''
		Initialize a Cipher instance with an alphabet of characters

		For messaeg-decoding purposes, best results found when
		alphabet is restricted to lowercase characters and space.
		This can be done with the clean_doc method
		'''
		self.alphabet = alphabet
		self.key=None

	def random_key(self):
		'''
		Generates a random trial key for the cipher
		'''
		#generates random trial key
		self.key={}
		unused=self.alphabet.copy()
		for c1 in self.alphabet:
			c2=random.choice(unused)
			self.key[c1]=c2
			unused.remove(c2)
		return None

	def transition_key(self,n=1):
		'''
		Transitions key randomly n (defaults to 1) number of times

		A transition involves two characters being chosen at random,
		and switching their corresponding cipher-mapped values.
		'''
		if self.key != None and len(self.key) > 2:
			for j in range(n):
				keys = list(self.key.keys())
				switched_keys = random.choices(keys,k=2)
				c1,c2 = switched_keys[0],switched_keys[1]
				self.key[c1], self.key[c2] = self.key[c2], self.key[c1]

	def encrypt(self, text):
		'''
		Encrypts input doc (a string) using the current key, returning encrypted string
		'''
		s=''
		for char in text:
			s+=self.key[char]
		return s

	def invert_key(self):
		'''
		Inverts the encryption key, converting the current key to its 
		corresponding solution key
		'''
		if self.key != None:
			pairs = [(c1, self.key[c1]) for c1 in self.key]
			for pair in pairs:
				self.key[pair[1]]=pair[0]

	def translate(self,text,read=True,key=None):
		'''
		Translates text using the current key, or a given key
		'''
		if key == None:
			key = self.key
		s = ''
		for j in text:
			s += key[j] 
		if read:
			print(s)
		return s

def log_likelihood(M,encrypted,cipher):
	'''
	Computes loglikelihood of a given key/cipher being tested

	Uses an (encrypted) string and corresponding test cipher (key) with 
	a transition frequency matrix M to evaluate logLikelihood for the cipher

	Parameters:
	M: first-order letter transition frequency matrix
	encrypted: a string which has been encrypted via substitution cipher
	cipher: a Cipher object used to translate the encrypted string into a 
	postulated message, for which the log likelihood (posterior prob.)
	is computed

	Returns:
	f: The logLikelihood of the cipher being tested, this is the log
	of the posterior probability of the input key given the observed
	encrypted string.
	'''
	f=0
	doc=cipher.translate(encrypted,read=False)
	doc_size=len(doc)
	for j in range(doc_size-1):
		c1, c2 = doc[j], doc[j+1]
		if c1 in M:
			if c2 in M[c1]:
				f+=np.log(M[c1][c2])
			else:
				f+=np.log(0.01*min([M[c1][k2] for k2 in M[c1]]))
	return f

def clean_doc(file):
	'''
	Converts csv file into string, removing punctuation and digits
	'''
	removed_chars=string.punctuation+string.digits
	translator = str.maketrans(removed_chars, len(removed_chars)*' ')
	s=''
	with open(file,'r') as input_text:
		reader=csv.reader(input_text)
		for row in reader:
			if len(row)>0:
				line=''.join(row)
				line=line.translate(translator).lower()
				s+=line+' '
	return s

def letter_freq_distribution(s_ref):
	'''
	Uses reference document s_ref to construct letter-transition frequency matrix

	Input:
	s_ref: a string free of punctuation and digits, intended to be long enough to
	allow accurate estimates of first order letter transitions

	Returns:
	M: a matrix (dictionary of dictionaries) whose entries are the 
	relative frequencies of a given letter transition, i.e. M['a']['b'']
	is the probability of observing a transition from 'a' to 'b''
	given that the first letter of the pair is 'a'
	'''
	M={}
	N={} #normalizing constants
	doc_size = len(s_ref)
	for i in range(doc_size-1):
		c1, c2 = s_ref[i], s_ref[i+1]
		if c1 in M:
			N[c1] += 1
			if c2 in M[c1]:
				M[c1][c2] += 1
			else:
				M[c1][c2] = 1
		else:
			M[c1]={}
			M[c1][c2] = 1
			N[c1] = 1
	for c1 in M:
		for c2 in M[c1]:
			M[c1][c2] = M[c1][c2]/N[c1]

	return M

def Decoder(M,encrypted,cipher,n,l=log_likelihood,progress=False,plot=False):
	'''
	Function uses MCMC simulation to solve substitution ciphers

	Uses the Metropolis Hastings algorithm to search the space of
	possible cipher solution keys, using the Cipher.transition method
	as the proposal function. The MAP estimator is returned as the 
	estimate of the decrypted text.

	Parameters:
	M: letter-transition frequency distribution
	encrypted: an encrypted (punctuation & digit)-free string
	cipher: a cipher object which evolves according to MH algorithm
	n: the maximum no. of accepted transtions. 1000 works well for medium-sized msgs
	l: log-likelihood function
	progress: if True, prints percentage to completion and current decryption state
	plot: if True, plots loglikelihood against iteration number

	Returns:
	decrypted: the most plausible estimate of the decrypted document found


	'''
	accepted_transitions = 0
	log_p = []
	l_max = l(M,encrypted,cipher)
	opt_key = cipher.key.copy()
	while accepted_transitions < n:
		k = cipher.key.copy()
		l_0 = l(M,encrypted,cipher)
		cipher.transition_key()
		l_trial = l(M,encrypted,cipher)
		dl = l_trial-l_0
		trans_prob = min(1,np.exp(dl))
		if random.random() < trans_prob:
			accepted_transitions += 1
		else:
			cipher.key = k		
		l_new = l(M,encrypted,cipher)
		if l_new > l_max:
			l_max = l_new
			opt_key = cipher.key.copy()
		log_p.append(l(M,encrypted,cipher))
		if accepted_transitions%10**1 == 0 and progress:
			print("~ ", 100*accepted_transitions/n, " pct completed.")
			cipher.translate(encrypted,key=opt_key)

	decrypted=cipher.translate(encrypted,key=opt_key,read=False)
	if plot:
		plt.plot([i+1 for i in range(len(log_p))],log_p)
		plt.xlabel("Iterations")
		plt.ylabel("Relative log-likelihood")
		plt.show()
	return decrypted

def init_guess(M, encrypted, cipher, n, m=3, l=log_likelihood):
	'''
	Computes a reliable initial guess key for the MH Decoder above

	Runs the decoder algorithm m times (default is 3)
	with max accepted iterations (n) intended to be considerably 
	smaller than that used in the actual Decoder.
	'''
	l_max = l(M, encrypted, cipher)
	opt_key = cipher.key
	for j in range(m):
		Decoder(M, encrypted, cipher, n)
		new_l=l(M, encrypted, cipher)
		if new_l>l_max:
			l_max, opt_key = new_l, cipher.key
		cipher.random_key()

	cipher.key=opt_key

if __name__ == '__main__':

	text_name = 'text_sample.csv'
	n=10**3

	reference_string = clean_doc('War-and-Peace.csv')
	M = letter_freq_distribution(reference_string)
	alphabet = list(M.keys())

	text_sample = clean_doc(text_name)

	myCipher=Cipher(alphabet)
	myCipher.random_key()
	encrypted=myCipher.encrypt(text_sample)
	myCipher.invert_key()
	solution=myCipher.key.copy()
	myCipher.random_key()

	print('Encrypted file is below: ')
	print(encrypted)
	print('--------------------------')
	decrypted=Decoder(M, encrypted, myCipher, n, progress=True)
	print("Here is the optimal text found: ")
	print('----------------------------')
	print(decrypted)