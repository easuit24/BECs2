import numpy as np

class Point(): 
    def __init__(self, xcoor, ycoor, vortex, trajectory = []): 
        self.xcoor = xcoor
        self.ycoor = ycoor
        self.vortex = vortex # true if vortex, false if anti-vortex
        self.trajectory = trajectory 

    def getCoors(self): 
        return (self.xcoor, self.ycoor) 
    
    def getVortexType(self): 
        return self.vortex 
    
    def getTrajectory(self): 
        return np.array(self.trajectory)  
    
    def addCoor(self, x, y): 
        self.xcoor = x
        self.ycoor = y
        self.trajectory.append((x,y))



class PointTracker(): 
    def __init__(self, psi_snaps, dx, L, points = [], border_threshold = 4): 
        self.points = points # an array of point objects 
        self.psi_snaps = psi_snaps
        self.dx = dx
        self.L = L 
        self.border_threshold = border_threshold
        self.circulation = None

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

        if len(avp) > 0: # initialize for anti-vortex 
            for i in range(len(avp)): 
                self.points.append(Point(avp[i][0], avp[i][1], trajectory = [(avp[i][0], avp[i][1])], vortex = False))

    def runTracker(self): 
        self.initGrid()

    def labelVortices(self): 
        print(len(self.psi_snaps))
        for i in range(1,len(self.psi_snaps)): 
            
            vortex_positions, anti_vortex_positions = self.detectVortices(self.psi_snaps[i])
            print("Outer Loop: ", i)
            for j in range(len(self.points)): 
                print("Inner Loop: ", i)
                # match the vortices in this array with the existing vortices in the self.points array 
                # start with the existing points and see which detected vortex is closest to each of the existing vortices

                # distance between detected points and existing point 
                # test this! I think this assumes that there are more detected points than existing points - the ADD point case 
                existing_point = self.points[j] 
                if existing_point.getVortexType() == True: # then it is a vortex 
                    detected_points = vortex_positions
                    euclidean_distances = np.abs(existing_point[0] - detected_points[:,0])**2 + np.abs(existing_point[1] - detected_points[:,1])**2
                    min_index = np.where(euclidean_distances == np.min(euclidean_distances))
                    min_coordinate = vortex_positions[min_index]
                else: 
                    detected_points = anti_vortex_positions
                    euclidean_distances = np.sqrt(np.abs(existing_point.getCoors()[0] - detected_points[:,0])**2 + np.abs(existing_point.getCoors()[1] - detected_points[:,1])**2)
                    min_index = np.where(euclidean_distances == np.min(euclidean_distances))
                    min_coordinate = anti_vortex_positions[min_index]
                existing_point.addCoor(*min_coordinate[0]) 
                self.points[j] = existing_point 
                

            if len(self.points) < len(vortex_positions) + len(anti_vortex_positions): 
                # initialize new points that correspond to the new vortices
                break  

            elif len(self.points) > len(vortex_positions) + len(anti_vortex_positions): 
                # end the points that correspond to the unclaimed points 
                break 
                 


