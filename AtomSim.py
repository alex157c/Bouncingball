import numpy as np
import pygame
import math
import os


class atoms:
    def __init__(self, x, y, ref, screen, mass, h, ground_y_top, k):
        self.obj_ref = ref
        #r is the position vector
        self.r = np.array([x, y])
        self.r_new = np.array([x, y])
        #Magnituides of the vectors to the neighbors at rest
        self.d = [None]*8
        #Reference to neighbors
        self.neighbors = [None]*8
        #Screen to draw on
        self.screen = screen
        #unit kg
        self.mass = mass
        #spring constant
        self.k = k
        #distance in pixels from the top of the window that the surface is located(same for ground_y_bottom just distance from bottom instead).
        #Thus needing to apply a force while under it. (only important for displaying since it wouldn't be visable outside the window)
        self.ground_y_top = ground_y_top
        #time step
        self.h = h
        #resulting force
        self.F_res = np.array([0.0, 0.0])
        #gravity force
        self.F_g = np.array([0.0, 0.0])
        #spring forces len=8 one for each neighbor
        self.F_s = np.array([[0.0, 0.0], 
                             [0.0, 0.0],
                             [0.0, 0.0],
                             [0.0, 0.0],
                             [0.0, 0.0],
                             [0.0, 0.0],
                             [0.0, 0.0],
                             [0.0, 0.0]])
        #spring force from the surface
        self.F_s_surface = np.array([0.0, 0.0])
        #acceleration
        self.a = np.array([0.0, 0.0])
        #velocity 
        self.v = np.array([0.0, 0.0])
    def print(self):
        return self.d, self.neighbors
    #where h is the time step
    def find_new_pos(self):
        #unit m*s^-2
        g = 9.82
        self.F_g[1] = self.mass*g
        self.F_res = self.F_g
        if (self.r[1] >= self.ground_y_top):
            x = self.ground_y_top-self.r[1]
            self.F_s_surface[1] = x*self.k
            self.F_res += self.F_s_surface
        for i in range(len(self.neighbors)):
            d = self.neighbors[i].Distance_from_point(self.r)
            delta_d = self.d[i]-d
            F_s_magnituide = (delta_d)*self.k
            self.F_s[i] = ((self.r-self.neighbors[i].GetPosition())/d)*F_s_magnituide
            self.F_res += np.array([self.F_s[i][0], self.F_s[i][1]])
        self.a = self.F_res/self.mass
        self.v += self.a*self.h
        self.r_new += self.v*self.h
    def update_pos(self):
        self.r = self.r_new
    def draw(self):
        pygame.draw.rect(self.screen, [0, 0, 0], [self.r[0], self.r[1], 1, 1])
    def distance_from_centroid(self, centroid):
        cr = centroid-self.r
        return cr
    def GetPosition(self):
        return self.r
    def Distance_from_point(self, x1):
        x2=self.r
        d=np.sqrt(((x2[0]-x1[0])**2)+((x2[1]-x1[1])**2))
        return d
    def FindNeighbour(self):
        Increment = 0
        #Defining Neighbours from the object itself
        NeighBourCoordinates = [[self.r[0]-1, self.r[1]-1], [self.r[0], self.r[1]-1], [self.r[0]+1, self.r[1]-1], 
                                [self.r[0]-1, self.r[1]], [self.r[0]+1, self.r[1]], 
                                [self.r[0]-1, self.r[1]+1], [self.r[0], self.r[1]+1], [self.r[0]+1, self.r[1]+1]]
        
        for i in range(len(self.obj_ref)):
            for j in range(len(NeighBourCoordinates)):
                bool_neighbor = self.obj_ref[i].GetPosition() == NeighBourCoordinates[j]
                if(bool_neighbor[0] == True and bool_neighbor[1] == True):
                    self.neighbors[Increment] = self.obj_ref[i]
                    self.d[Increment] = self.obj_ref[i].Distance_from_point(self.r)
                    Increment += 1
        b = np.argwhere(self.neighbors)
        self.neighbors = self.neighbors[np.min(b):np.max(b)+1]
    def temp_draw(self):
        pygame.draw.rect(self.screen, [0, 255, 0], [self.r[0], self.r[1], 1, 1])


#n: grid of n*n points
#mass: mass of each point [kg]
#fps: frames per second
#runtime: amount of frames to run for
def setup(k=10, n=100, screen_size=(600, 600), mass=100, fps=60, runtime=100, ground_y_bottom=50):
    #pre-calculated values
    #initialize pygame
    pygame.init()
    #create a window
    screen = pygame.display.set_mode(screen_size)
    pygame.display.set_caption("Atom Sim")
    #clock is used to set a max fps
    clock = pygame.time.Clock()
    #radius of the ball of points
    r=(n/2)-1
    #creates a list of n^2 that will contain references to the atom objects
    obj_ref = [None]*(n*n)
    #time step
    #not sure this is correct (times 8 is just to make it quicker)
    h = (1/fps)*8
    #creates all the points
    #2xfor is slow. Try to figure a better method. Especially since it will be 3xfor in 3d. 
    #Need all points in a square with a homogeneous spacing
    for i in range(n):
        for j in range(n):
            #(screen_size[1]/2) makes it start the middle of the screen
            obj_ref[i*n+j] = atoms((screen_size[0]/2)-r+i, (screen_size[1]/2)-r+j, obj_ref, screen, mass, h, (screen_size[1]-ground_y_bottom), k)
    #removing points not in circle
    start_centroid = np.array([screen_size[0]/2, screen_size[1]/2])
    remove_indicies_obj_ref = [None]*len(obj_ref)
    for i in range(len(obj_ref)):
        dis_centroid = obj_ref[i].distance_from_centroid(start_centroid)
        if (math.sqrt(pow(dis_centroid[0], 2)+pow(dis_centroid[1], 2)) > r):
            remove_indicies_obj_ref[i] = i
    remove_indicies_obj_ref = list(filter(lambda num: num != None, remove_indicies_obj_ref))
    remove_indicies_obj_ref.sort(reverse=True)
    for i in range(len(remove_indicies_obj_ref)):
        del obj_ref[remove_indicies_obj_ref[i]]
    
    for i in range(len(obj_ref)):
        obj_ref[i].FindNeighbour() 
    npy_names = os.listdir("D:/Python Projects/Test Projects/AtomSim/Billeder")
    for i in range(len(npy_names)):
        os.remove(f'D:/Python Projects/Test Projects/AtomSim/Billeder/{npy_names[i]}')
    #saving loop
    system_state = [None]*len(obj_ref)
    for j in range(runtime):
        for i in range(len(obj_ref)):
            obj_ref[i].find_new_pos()
        for i in range(len(obj_ref)):
            obj_ref[i].update_pos()
            system_state[i] = obj_ref[i].GetPosition()
        np.save(f'D:/Python Projects/Test Projects/AtomSim/Billeder/{j}', system_state, allow_pickle=True)
    draw(screen, clock, fps)


def draw(screen, clock, fps):
    npy_names = os.listdir("D:/Python Projects/Test Projects/AtomSim/Billeder")
    max_frames = len(npy_names)
    running = True
    incrementor = 0
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
        #clear the screen
        screen.fill([255, 255, 255])
        system_state = np.load(f'D:/Python Projects/Test Projects/AtomSim/Billeder/{incrementor}.npy', allow_pickle=True)
        #main loop for drawing
        for i in range(len(system_state)):
            pygame.draw.rect(screen, [0, 0, 0], [*system_state[i], 1, 1])
        #pygame.draw.rect(screen, [0, 0, 0], [350, 350, 1, 1])
        pygame.display.flip()
        clock.tick(fps)
        if (incrementor != max_frames-1):
            incrementor += 1
        else:
            incrementor = 0

#saving and displaying
setup(n=80, screen_size=(800, 800))