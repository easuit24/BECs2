# New tracking algorithm for tracking a pair of vortices 

import numpy as np
from scipy.cluster.vq import vq 
from collections import defaultdict 

class Point(): 
    def __init__(self, xcoor, ycoor, vortex, trajectory = [], starttime = 0): 
        self.xcoor = xcoor
        self.ycoor = ycoor
        self.vortex = vortex # true if vortex, false if anti-vortex
        self.trajectory = trajectory 
        self.starttime = starttime 
        self.endtime = None 

    def getCoors(self): 
        return (self.xcoor, self.ycoor) 
    
    def getVortexType(self): 
        return self.vortex 
    
    def getTrajectory(self): 
        return np.array(self.trajectory)  
    
    def getStartTime(self): 
        return self.starttime
    
    def calcIndexBuffer(self): 
        '''
        The starting index where the vortex is introduced 
        '''
        return self.starttime/(250*self.dt) 
        
    
    def addCoor(self, x, y): 
        self.xcoor = x
        self.ycoor = y
        self.trajectory.append((x,y))

    def endPoint(self, x, y, end): 
        self.endtime = end
        self.xcoor = x
        self.ycoor = y 


class PointTracker(): 
    def __init__(self, psi_snaps, dx, L, dt, points = [], border_threshold = 4): 
        self.points = points # an array of currently active point objects 
        self.vortices = [] 
        self.antivortices = [] 
        self.point_history = points # all points that were active during the simulation
        self.unorderedpoints = points 
        self.psi_snaps = psi_snaps
        self.dx = dx
        self.L = L 
        self.dt = dt 
        self.border_threshold = border_threshold
        self.circulation = None
        self.trajectories = []

        self.runTracker()
    
    def getPoints(self):
        return self.points 
    
    def detectVortices(self, psi):
        # extract the inner part of the box
        psi = psi[int(self.L/2/self.dx):int(3*self.L/2/self.dx), int(self.L/2/self.dx):int(3*self.L/2/self.dx)]
    
        Nx, Ny = np.shape(psi)
        S = np.angle(psi)
        vortex_positions = [] 
        anti_vortex_positions = [] 
    
        # initialize arrays
        dS_y_left = np.zeros((len(psi[0]), len(psi[0])))
        dS_y_right = np.zeros((len(psi[0]), len(psi[0])))
        dS_x_top = np.zeros((len(psi[0]), len(psi[0])))
        dS_x_bottom = np.zeros((len(psi[0]), len(psi[0])))
    
        for i in range(self.border_threshold, Nx-self.border_threshold):
            for j in range(self.border_threshold,Ny-self.border_threshold):
                dS_y_left[i,j] = np.mod((S[i, j+1]-S[i,j])+np.pi, 2*np.pi)-np.pi
                dS_y_right[i,j] = np.mod((S[i+1, j+1] - S[i+1, j])+np.pi, 2*np.pi)-np.pi
                dS_x_top[i,j] = np.mod((S[i+1, j+1]-S[i,j+1])+np.pi, 2*np.pi)-np.pi
                dS_x_bottom[i,j] = np.mod((S[i+1,j]-S[i,j])+np.pi, 2*np.pi)-np.pi
                circulation_ij = -dS_y_left[i,j] -dS_x_top[i,j] + dS_y_right[i,j] + dS_x_bottom[i,j] 

                if circulation_ij > 6.2: 
                    vortex_positions.append([(j+0.5)*self.dx, (i+0.5)*self.dx])
                elif circulation_ij < -6.2: 
                    anti_vortex_positions.append([(j+0.5)*self.dx, (i+0.5)*self.dx]) 

        circulation  = -dS_y_left -dS_x_top + dS_y_right + dS_x_bottom  
        

        # plot the circulation to see where the vortices are found!
        self.circulation = circulation 
        return np.array(vortex_positions), np.array(anti_vortex_positions)
    
    def initGrid(self): 
        vp, avp = self.detectVortices(self.psi_snaps[0])
 
        if len(vp) > 0: # initialize for vortex
            for i in range(len(vp)): 
                self.points.append(Point(vp[i][0], vp[i][1], trajectory = [(vp[i][0], vp[i][1])], vortex = True)) 
                self.vortices.append(Point(vp[i][0], vp[i][1], trajectory = [(vp[i][0], vp[i][1])], vortex = True)) 

        if len(avp) > 0: # initialize for anti-vortex 
            for i in range(len(avp)): 
                self.points.append(Point(avp[i][0], avp[i][1], trajectory = [(avp[i][0], avp[i][1])], vortex = False))
                self.antivortices.append(Point(avp[i][0], avp[i][1], trajectory = [(avp[i][0], avp[i][1])], vortex = False))

    def getCurrentPointInfo(self): 
        self.trajectories = [] 
        #self.vtraj = [] 
        #self.avtraj = [] 
        self.vtypes = [] 
        for i in range(len(self.points)):

            self.trajectories.append(self.points[i].getCoors())
            #self.vtypes.append(self.points[i].getVortexType())
            #self.vtraj.append(self.vortices[i].getCoors())
            self.avtraj.append(self.antivorticesvortices[i].getCoors())
        return self.trajectories, self.vtypes 
    
    def getVAntiVInfo(self): 
        self.vtraj = [] 
        self.avtraj = []
        for i in range(len(self.vortices)): 
            self.vtraj.append(self.vortices[i].getCoors())
        for i in range(len(self.antivortices)): 
            self.avtraj.append(self.antivortices[i].getCoors())

        return self.vtraj, self.avtraj
        


    def runTracker(self): 
        self.initGrid()

    def labelVortices(self): 

        # loop through each of the frames of the simulation 
        for i in range(1,len(self.psi_snaps)):  
            # find the vortex and anti-vortex positions           
            vortex_positions, anti_vortex_positions = self.detectVortices(self.psi_snaps[i]) # this can be the features
            
            # find active point coordinates 
            active_vortex_coors, active_antivortex_coors = self.getVAntiVInfo() # this can be the codebook 

            vortices_closest_index, distances_v = vq(active_vortex_coors, vortex_positions)
            antivortices_closest_index, distances_av = vq(active_antivortex_coors, anti_vortex_positions)

            print("Active Points: ", active_vortex_coors + active_antivortex_coors)
            print("New Detected Vortices: ", vortex_positions)
            print("New Detected AntiVortices: ", anti_vortex_positions)
            print('Closest Match: ' , vortices_closest_index)
            print('AntiVortex Closest Match: ', antivortices_closest_index)





                        