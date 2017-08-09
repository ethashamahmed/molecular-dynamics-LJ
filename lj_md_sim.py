import numpy as np 

def md():
    # main function of the molecular dynamics program

    global file_xyz
    
    initialisation().variables()
    file_xyz = open('output.xyz', 'w')

    print ' '
    print 'This is a molecular dynamics program'
    print ' '
    print 'Step       Temperature      Potential       Kinetic        Total '
    print '           (K)              Energy (J)      Energy (J)     Energy(K) ' 

    for step in range(Nsteps+1):
        if step == 0:
            initialisation().initial_pos(cluster1)
            initialisation().initial_pos(cluster2)
            initialisation().center_of_mass()
            initialisation().inital_velocity()
            compute().lennard_jones()
            compute().temperature()
            output(step)
        
        else:
            compute().temperature()

            if step % 10 ==0:
                #prints a frame every 10 steps to .xyz file 
                output(step) 

            update().positions()
            update().velocities()
            compute().lennard_jones()
            update().velocities()
            compute().total_energies()

            print step, temp, ene_pot_tot, ene_kin_tot, ene_tot

    file_xyz.close()

class initialisation():
    def variables(self):

        global Nsteps, dt, Rcutoff, phicutoff, init_distance, init_speed, impact_parameter, Ntotal, DIM, mass, sigma, ep, kB,pos, vel, acc, Cmass, sqspeed, ene_pot, ene_kin, c1, c2, cluster1, cluster2, particle, conver
        
        #initial parameters
        Nsteps = 30000              # Nsteps needs to be a very large integer because dt by default is small. 
        dt = 1.0e-4                 # Time step, must be small enough to calculate accurately but large enough to run similation at reasonable speed. 
        conver = 1.0e-12/dt         # used to convert magnitude of temeperature to kelvin (magnitude of dt and conver must add up to 1e-12)
        init_speed = 2.0            # initial speed of the system. 2.0 by default 
        init_distance = 5.0         # initial distance between the clusters. 5.0 by default 
        impact_parameter = 1.0      # shifts cluster two in y-direction 
        Rcutoff = 2.5               # cut off distance for lenard jones potential 
        phicutoff = 4 * ((Rcutoff**-12) - (Rcutoff ** -6))  #used to keep potential at Rcutoff 0

        #constants used to convert to SI 
        #mass, sigma, and ep have been chosen for Ar
        mass = 6.634e-26            # mass in kg
        sigma = 3.4e-10             # equilibrium distance between two particles in m 
        ep = 1.653e-21              # depth of potential well when r = sigma in joules
        kB = 1.381e-23              # Boltzmann Constant in m^2 kg s^-2 K^-1
               

        #particle parameters
        c1 = int(input('Choose number of particles in cluster 1: '))        # size of cluster 1 from user input
        c2 = int(input('Choose number of particles in cluster 2: '))        # size of cluster 2 from user input 
        particle = 'Ar'            # name of particle used
        Ntotal = c1 + c2           # total number of particles
        DIM = 3                    # Spatial dimensions 

        #initialising arrays
        pos = [[], [], []]                                                  # position array
        vel = [np.zeros(Ntotal),np.zeros(Ntotal),np.zeros(Ntotal)]          # velocity array 
        acc = [[], [], []]                                                  # acceleration array 
        Cmass = np.zeros(DIM)                                               # center of mass array 
        sqspeed = np.zeros(Ntotal)                                          # velocity squared array 
        ene_pot = []                                                        # potential energy array 
        ene_kin = np.zeros(Ntotal)                                          # kinetic energy array 

        #input filename allocation  
        cluster1 = raw_input('Choose input filename for cluster 1: ')
        cluster2 = raw_input('Choose input filename for cluster 2: ')

    def initial_pos(self,cluster):
        # reads initial coordinates from file
        global pos

        file = open(cluster, 'r')
        for line in file.readlines():
            xyz = line.split()
            for k in range(DIM):
                # Inserts the initial coordinates into pos array
                pos[k].append(float(xyz[k]))
        file.close()
    
    def center_of_mass(self):
        #calculates center of mass of the system 
        global pos

        for k in range(DIM):
                Cmass[k] = np.mean(pos[k][:c1])     #calculates center of mass of cluster 1
        init_pos = -init_distance*(float(c2)/float(Ntotal))

        #translates cluster 1
        pos[0][:c1] = [(value - Cmass[0] + init_pos) for value in pos[0][:c1]]
        pos[1][:c1] = [(value - Cmass[1]) for value in pos[1][:c1]]
        pos[2][:c1] = [(value - Cmass[2]) for value in pos[2][:c1]]

        for k in range(DIM):
            Cmass[k] = np.mean(pos[k][c1:])         #calculates center of mass of cluster 2
        init_pos = init_distance*(float(c1)/float(Ntotal))

        #translates cluster 2
        pos[0][c1:] = [(value - Cmass[0] + init_pos) for value in pos[0][c1:]]
        pos[1][c1:] = [(value - Cmass[1] + float(impact_parameter)) for value in pos[1][c1:]]
        pos[2][c1:] = [(value - Cmass[2]) for value in pos[2][c1:]]


    def inital_velocity(self):
        #velocity initialised
        global vel

        #keeps center of mass constant
        init_vel_1 = init_speed*(float(c1)/float(Ntotal))
        #opposite direction
        init_vel_2 = (-init_speed)*(float(c2)/float(Ntotal))

        #updates initial velocity in x direction 
        vel[0][:c1] = [(value + init_vel_1) for value in vel[0][:c1]]
        vel[0][c1:] = [(value + init_vel_2) for value in vel[0][c1:]]

class compute():
    def temperature(self):
        # Temperature is calculated in Kelvin 
        global ene_kin, temp

        for i in range(Ntotal):
            sqspeed[i] = vel[0][i]**2 + vel[1][i]**2 + vel[2][i]**2
            ene_kin[i] = sqspeed[i]*0.5     #kinetic energy calculated and stored in ene_kin array
        mean_sq_speed = np.mean(sqspeed)    #mean square speed
        temp = ((mean_sq_speed)/(DIM*kB))* mass * (sigma/(dt*conver))**2      #temperature calculated using equipartition theorem and (sigma/dt*conver)**2 is used to convert reduced unit back to SI

    def lennard_jones(self):
        #calculated the forces due to Lennard Jones potential
        global ene_pot, acc

        #initialises and resets force, potential and acceleration 
        r = np.zeros(DIM)
        force = np.zeros(DIM)
        ene_pot = np.zeros(Ntotal)
        acc = [np.zeros(Ntotal),np.zeros(Ntotal),np.zeros(Ntotal)]        

        #looping over every pair of particles
        for i in range(Ntotal):
            for j in range(Ntotal):
                if (i <= j):        #excludes paris that have already been calculated and the same particles
                    continue
                else:
                    for k in range(DIM):        #particle pair seperation
                        r[k] = pos[k][i] - pos[k][j]

                    r_abs = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
                    
                    if r_abs < Rcutoff:     #checks if the particles are within the cutoff distance

                        #calculated Lennard Jones potential
                        lj_pot = ( (4*( (r_abs**-12) - (r_abs**-6) ) ) - phicutoff)
                        #updates potential energy array
                        ene_pot[i] += lj_pot
                        ene_pot[j] += lj_pot

                        #F = -div(V)
                        lj_force = 24 * ((2*r_abs**-13)- (r_abs**-7))

                        for k in range(DIM):
                            #Fx = dV/dr * dr/dx
                            #dr/dx = x/r 
                            force[k] = lj_force * (r[k]/r_abs)

                            #accelaration array updated using force due to Lennard jones potential
                            #Fij = -Fji
                            # a(t+dt) = f(t) / m where m = 1 in reduced units 
                            acc[k][i] += force[k]
                            acc[k][j] -= force[k]

    def total_energies(self):
        global ene_kin_tot, ene_pot_tot, ene_tot 
        #total kinetic, total potential and total energies calculated and converted to SI from reduced unit
        ene_kin_tot = sum(ene_kin) * mass 
        ene_pot_tot = sum(ene_pot) * ep
        ene_tot = ene_kin_tot + ene_pot_tot

class update():
    #velocity verlet algorithm 

    def positions(self):
        # x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
        global vel, acc
        for i in range(Ntotal):
            for k in range(DIM):
                pos[k][i] += vel[k][i]*dt + 0.5 * acc[k][i] * dt**2
        
    def velocities(self):
        # v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
        global vel
        for i in range(Ntotal):
            for k in range(DIM):
                vel[k][i] += 0.5 * acc[k][i] * dt

def output(step):
    #writes to a .xyz file to be used in a visualisation program like VMD 
    Time = dt*step  #calculates time
    file_xyz.write(str(Ntotal) +'\n')
    file_xyz.write('Step: %s	Time: %sps	Temperature: %sK\n' % (step, Time, temp))
    for i in range(Ntotal):
		file_xyz.write('%s	%s	%s	%s\n' % (particle, pos[0][i], pos[1][i], pos[2][i]))



md()
