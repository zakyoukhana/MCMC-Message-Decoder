import csv
import random
import numpy as np
import os
import matplotlib.pyplot as plt
import string

csv.register_dialect('myDialect',
delimiter = ',',
skipinitialspace=True)

class cipher(object):
	#substitution cipher object
	def __init__(self, alphabet):
		self.alphabet=list(set(alphabet))
		self.key=None

	def random_key(self):
		#generates random trial key
		self.key={}
		unused=self.alphabet.copy()
		for c1 in self.alphabet:
			c2=random.choice(unused)
			self.key[c1]=c2
			unused.remove(c2)
		return None

	def transition_key(self,n=1):
		#transitions key n number of times
		k_v=False
		if self.key!=None and len(self.key)>2:
			for j in range(n):
				keys=list(self.key.keys())
				switched_keys=random.choices(keys,k=2)
				c1,c2=switched_keys[0],switched_keys[1]
				self.key[c1],self.key[c2]=self.key[c2],self.key[c1]
		return None

	def encrypt(self, file):
		#encrypts file using current key
		with open('encrypted_file.csv', mode='w') as encrypted, open(file, 'r') as csvFile:
			encrypter = csv.writer(encrypted)
			reader = csv.reader(csvFile, dialect='myDialect')
			for row in reader:
				if len(row)>0:
					phrase=row[0]
					encrypter.writerow([''.join([self.key[c] for c in phrase]),])
		return 'encrypted_file.csv'

	def invert_key(self):
		#inverts key, returning solution for the current key
		if self.key!=None:
			pairs=[(c1,self.key[c1]) for c1 in self.key]
			for pair in pairs:
				self.key[pair[1]]=pair[0]
		return None

	def translate(self,text,read=True,key=None):
		#translates an encrypted text using current key (or a given key)
		if key==None:
			key=self.key
		with open('decrypted_file.csv', mode='w') as decrypted, open(text, 'r') as csvFile:
			decrypter = csv.writer(decrypted)
			reader = csv.reader(csvFile, dialect='myDialect')
			for row in reader:
				if len(row)>0:
					phrase=row[0]
					n=len(phrase)
					decrypter.writerow([''.join([key[c] for c in phrase]),])
		if read:
			with open('decrypted_file.csv','r') as decrypted:
				reader=csv.reader(decrypted)
				for row in reader:
					print(row)
		return 'decrypted_file.csv'

	def init_guess(self,M,l,encrypted):
		#comes up with initial guess key for MH algo thru large transitions
		if self.key!=None:
			l_max=l(M,encrypted,self.key)
			max_key=self.key.copy()
			for i in range(10):
				self.transition_key(n=10)
				l_trial=l(M,encrypted,self.key)
				if l_trial>l_max:
					l_max=l_trial
					max_key=self.key.copy()
			self.key=max_key
		return None

def logLikelihood(M,text_file,key):
	#computes loglikelihood of a given key
	f=0
	with open(text_file, 'r') as csvFile:
	    reader = csv.reader(csvFile, dialect='myDialect')
	    for row in reader:
	    	if len(row)>0:
	    		phrase=row[0]
	    		num_chars=len(phrase)
	    		for i in range(num_chars-1):
	    			c1,c2=key[phrase[i]],key[phrase[i+1]]
	    			if c1 in M:
	    				if c2 in M[c1]:
	    					f+=np.log(M[c1][c2])
	    				else:
	    					f+=np.log(0.01*min([M[c1][k2] for k2 in M[c1]]))
	return f

def clean_doc(file):
	translator = str.maketrans('', '', string.punctuation+string.digits)
	with open(file,'r') as input_text,open('clean.csv','w') as clean:
		reader=csv.reader(input_text)
		writer=csv.writer(clean)
		for row in reader:
			if len(row)>0:
				line=[row[0].translate(translator).lower()]
				writer.writerow(line)
	os.remove(file)
	os.rename('clean.csv',file)
	return None


alpha_1=set([]) #alphabet of all characters from cleaned war and peace
M={} #transition probability matrix
N_c1={}
#create transition frequency matrix
clean_doc('War-and-Peace.csv')
with open('War-and-Peace.csv', 'r') as csvFile:
    reader = csv.reader(csvFile, dialect='myDialect')
    for row in reader:
    	if len(row)>0:
    		phrase=row[0]
    		num_chars=len(phrase)
    		for i in range(num_chars-1):
    			c1, c2 =phrase[i], phrase[i+1]
    			if c1 in M:
    				N_c1[c1]+=1
    				if c2 in M[c1]:
    					M[c1][c2]+=1
    				else:
    					M[c1][c2]=1
    			else:
    				M[c1]={}
    				M[c1][c2]=1
    				alpha_1.add(c1)
    				N_c1[c1]=1
for c1 in M:
	for c2 in M[c1]:
		M[c1][c2]=M[c1][c2]/N_c1[c1]

text_name='text_sample.csv'
clean_doc(text_name)

#War&Peace excerpt to encrypt and test on
# with open(text_name,'w') as text,open('War-and-Peace.csv','r') as book:
# 	writer=csv.writer(text)
# 	reader=csv.reader(book)
# 	for i, row in enumerate(reader):
# 		if i>850 and len(row)>0:
# 			phrase=row[0]
# 			print(phrase)
# 			writer.writerow([phrase,])
# 		if i>875:
# 			break

myCipher=cipher(alpha_1)
myCipher.random_key()
encrypted=myCipher.encrypt(text_name)
myCipher.invert_key()
solution=myCipher.key.copy()
myCipher.random_key()

def Decoder(M,encrypted=encrypted,cipher=myCipher,l=logLikelihood,n=6*10**2):
	#MH Decoding algorithm implementation
	accepted_transitions=0
	myCipher.init_guess(M,l,encrypted)
	print(l(M,encrypted,cipher.key)," initial guess logLikelihood")
	log_p=[]
	l_max=l(M,encrypted,cipher.key)
	opt_key=cipher.key.copy()
	while accepted_transitions<n:
		k=cipher.key.copy()
		l_0=l(M,encrypted,k)
		cipher.transition_key()
		trial_k=cipher.key.copy()
		l_trial=l(M,encrypted,trial_k)
		dl=l_trial-l_0
		trans_prob=min(1,np.exp(dl))
		if random.random()<trans_prob:
			accepted_transitions+=1
		else:
			cipher.key=k		
		l_new=l(M,encrypted,cipher.key)
		if l_new>l_max:
			l_max=l_new
			opt_key=cipher.key.copy()
		log_p.append(l(M,encrypted,cipher.key))
		if accepted_transitions%10**1==0:
			print("~ ", 100*accepted_transitions/n, " pct completed.")
			cipher.translate(encrypted,key=opt_key)

	# print('best key tried i.e. MAP esimator')
	decrypted=cipher.translate(encrypted,key=opt_key,read=False)

	l_sol=l(M,encrypted,solution)
	plt.plot([i+1 for i in range(len(log_p))],log_p)
	plt.xlabel("Iterations")
	plt.ylabel("Relative log-likelihood")
	plt.show()
	return decrypted

decrypted=Decoder(M)
print("Here is the optimal text found: ")
print(open(decrypted).read())

#note: possibility to opt performance, adjust hyper-parameter dynamically