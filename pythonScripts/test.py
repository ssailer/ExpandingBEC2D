
#!/usr/bin/env python
import numpy as np
import math
# import h5py
# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.cluster.vq import kmeans, kmeans2, whiten


def main():



	threads = 16
	frontx = [None] * threads
	endx = [None] * threads
	partx = 4096 / threads
	for i in range(0,threads):
		if(i == 0):
			frontx[i] = (i * partx) +1
		else:
			frontx[i] = (i * partx)

		if((i == threads-1) | (i == 0)):
			endx[i] = partx-1
		else:
			endx[i] = partx
	index1 = frontx
	index2 = []
	for i in range(0,len(frontx)):
		index2.append(frontx[i] + endx[i] -1)

	print index1
	print index2

	print partx
	print frontx
	print endx
	# for(int i = 0; i < threads; i++){
		# if(i == 0){ frontx[i] = (i * partx) + 1;}
		# else{ frontx[i] = (i *partx);}
		# if(i == threads-1){ endx[i] = partx-1;}
		# else{endx[i] = partx;}
	# }
	 # for (int i = 0; i < threads; ++i){
	 	# cerr << "Thread# " << omp_get_thread_num() << endl;
	 	# cerr << "kblock = " << frontx[i] << "," << 1 << "," << endx[i] << "," << suby << endl;
	 	# cerr << frontx[i]-1 << "," << 1 << "," << endx[i] << "," << suby << endl;
	 	# cerr << frontx[i]   << "," << 1 << "," << endx[i] << "," << suby << endl;
	 	# cerr << frontx[i]+1 << "," << 1 << "," << endx[i] << "," << suby << endl;
	 	# cerr << frontx[i]   << "," << 0 << "," << endx[i] << "," << suby << endl;
	 	# cerr << frontx[i]   << "," << 1 << "," << endx[i] << "," << suby << endl;
	 	# cerr << frontx[i]   << "," << 2 << "," << endx[i] << "," << suby << endl;
	 # }


if __name__ == '__main__':
   	main()

