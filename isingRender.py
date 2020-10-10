#Created by Murat Magomedov on 7/29/19

#Simply place in the same directory as the '/output' folder and run 

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.colors as c
import imageio
import os
from tqdm import tqdm       #creates fancy loading bar if wrapped around an iterable object

def format(txtFile):        
    #creates a 2D list of integer colors based on text data for a single Ising iteration
    mesh = np.zeros((104, 104))
    for y, row in enumerate(txtFile):
        #cleans up text data
        row = row.replace('-1', '-')
        row = row.replace(' ', '') 

        for x, char in enumerate(row):
            if char == '1':
                mesh[y][x] = 1
        y += 1
    return mesh

def draw(mesh, imagePath):
	#matplotlib stuff to create an image based on the colormap for a single Ising iteration
    plt.figure(figsize=(5, 5), dpi=300)
    plt.axis('off')                         
    plt.pcolormesh(mesh, cmap = plt.cm.Greys)                            #creates the visual representation; modify 'cmap' to use other color schemes 
    plt.savefig(imagePath, bbox_inches='tight', pad_inches = 0)
    plt.close()                                                         #clears figure from memory after it has been saved 

def renderAll():
	#must have an /output/ folder containing text files in the same directory as this python file 
	#creates images in /states/ corresponding to every file in /output/
    print('Processing images...')
    if not os.path.exists('./states'):
        os.mkdir('./states')
    for filename in tqdm(sorted(os.listdir('./output'))):           #iterates through text files in /output/ directory; tqdm simply creates the loading bar in terminal
        file = open('./output/'+filename, 'r')                      #opens .txt in read only mode 
        colormap = format(file)             
        draw(colormap, './states/' + os.path.splitext(filename)[0]+'.png')
        file.close()
		
def gif():                                                              
    #creates a gif from images located in /states/ directory
	#pretty hackish and will fail if images aren't exactly named 'testX.png' or 'end_test.png'
    print('\nCreating gif...')
    images = []
    numImages = len(os.listdir('./states'))
    for i in tqdm(range(numImages-1)):                                                     
        images.append(imageio.imread('./states/test' + str(i) + '.png'))
    images.append(imageio.imread('./states/end_test.png'))
    imageio.mimsave('./ising.gif', images)

def mean(stableState):
    #visualizes arithmetic mean of states after stableState number of iterations 
    print("Averaging states...")
    numFiles = len(os.listdir('./output')) - 1
    meanCMap = np.zeros((104, 104))
    for i in tqdm(range(stableState, numFiles)):
        file = open('./output/test' + str(i) + '.txt')
        meanCMap += format(file)
    meanCMap /= (numFiles - stableState)
    for k in range(0, 2):
        meanCMap = np.delete(meanCMap, [0, 1, 102, 103], k)
    print(meanCMap.shape)
    draw(meanCMap, 'meanState3')

renderAll()
gif()
mean(500)