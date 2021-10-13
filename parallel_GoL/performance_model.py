import numpy as np
import matplotlib.pyplot as plt


x_grid = []
y_grid = []


figureHeight=3
figureWidth=3

plt.figure(figsize=(figureWidth,figureHeight))

panelWidth=5
panelHeight=3

relativePanelWidth=panelWidth/figureWidth
relativePanelHeight=panelHeight/figureHeight

panel=plt.axes([0.1,0.1,relativePanelWidth,relativePanelHeight])

panel.tick_params(bottom=True, labelbottom=True,
                   left=True, labelleft=True,
                   right=False, labelright=False,
                   top=False, labeltop=False)
panel.set_xlim(0,10)
panel.set_ylim(0,1.6)
panel.set_xlabel('Number of Processors')
panel.set_ylabel('Time (seconds)')
#panel.set_xticks([-10,-5,0,5,10])
#panel.set_yticks([0,10,20,30,40,50,60])



with open('data.txt', 'r') as infile:
    for line in infile:
        splitLine = line.strip().split('\t')
        #print (splitLine)
        aX = int(splitLine[0])
        aY = float(splitLine[1])
        x_grid.append(aX)
        y_grid.append(aY)


panel.plot(x_grid,y_grid, marker ='o', linewidth=0, markersize = 1)
p = np.polyfit(x_grid, y_grid, 3)


xLine = np.linspace(1,8,100)
yLine = np.linspace(1,8,100)


for elt in range(0,len(xLine)):
    yLine[elt] = p[0]*(xLine[elt]**3) + p[1]*(xLine[elt]**2) + p[2]*(xLine[elt]) + p[3]


panel.plot(xLine,yLine, 'r', linewidth=0.9)

t_c = 0.00227166132
t_s = 1.2867331327706748e-05
t_w = 2.5579058611622674e-09

xLine = np.linspace(1,8,100)
yLine = []

def T_func(P):
    T = (1/P)*(2*24*t_w + t_c*(24)**2) + 2*t_w*26 + 4*t_s
    return T

    
for P in xLine:
    yLine.append(T_func(P))

panel.plot(xLine,yLine, 'b', linewidth=0.9)
