import math
import sys
import time

starttime = time.time()

status = 0

def guage():
    sys.stdout.write('\r')
    sys.stdout.write('['+'=' * int(status/5)+' '*(20-int(status/5))+']{}%'.format(status))
    sys.stdout.flush()

def tuple_input(text):
    first_devision = text.split('(')[1]
    second_devision = first_devision.split(')')[0]
    third_devision = second_devision.split(',')
    to_return = (float(third_devision[0]), float(third_devision[1]), float(third_devision[2]))
    return to_return

def distance_2D(a, b):
    x1 = a[0]
    y1 = a[1]
    x2 = b[0]
    y2 = b[1]

    X = x1-x2
    Y = y1-y2
    X2 = X**2
    Y2 = Y**2
    distance2 = X2+Y2
    distance = distance2**0.5

    return distance

g = 9.81

sleep_time = 0.01
motors_first_deg_unit = 90
motors_first_deg = []
for i in range(12):
    motors_first_deg.append(motors_first_deg_unit)

old_deg = motors_first_deg

motors_dir = [1, 1, 1, 1, 1, 1, 1, -1, -1, -1, 1, 1]

one_action_time = float(input('One Action Time:'))

body_mass = float(input('Body Particle Mass:'))

leg1_upper_particle = tuple_input(input('Leg1 Upper Particle:'))
leg1_knee_particle = tuple_input(input('Leg1 Knee Particle:'))
leg1_ancle_particle = tuple_input(input('Leg1 Ancle Particle:'))
leg2_upper_particle = tuple_input(input('Leg2 Upper Particle:'))
leg2_knee_particle = tuple_input(input('Leg2 Knee Particle:'))
leg2_ancle_particle = tuple_input(input('Leg2 Ancle Particle:'))
body_weight_particle = tuple_input(input('Body Particle:'))

below_leg = distance_2D((leg1_ancle_particle[2], leg1_ancle_particle[1]), (leg1_knee_particle[2], leg1_knee_particle[1]))
upper_leg = distance_2D((leg1_knee_particle[2], leg1_knee_particle[1]), (leg1_upper_particle[2], leg1_upper_particle[1]))

leg1_upper_mass = float(input('Leg1 Upper Mass:'))
leg1_knee_mass = float(input('Leg1 Knee Mass:'))
leg1_ancle_mass = float(input('Leg1 Ancle Mass:'))
leg2_upper_mass = float(input('Leg2 Upper Mass:'))
leg2_knee_mass = float(input('Leg2 Knee Mass:'))
leg2_ancle_mass = float(input('Leg2 Ancle Mass:'))

forwarding_rate = float(input('Forwarding Rate:'))
gap_allow = float(input('Gap allow:'))

file_dir = input('File Save Dir:')
file = open(file_dir, 'a')

print('\n====================')
print('Calculation Process\n====================\n')
guage()

file.write('leg1 particles:{}\n'.format([leg1_upper_particle, leg1_knee_particle, leg1_ancle_particle]))
file.write('leg2 particles:{}\n'.format([leg2_upper_particle, leg2_knee_particle, leg2_ancle_particle]))
file.write('central particle:{}\n'.format(body_weight_particle))
file.write('leg1 masses:{}\n'.format([leg1_upper_mass, leg1_knee_mass, leg1_ancle_mass]))
file.write('leg2 masses:{}\n'.format([leg2_upper_mass, leg2_knee_mass, leg2_ancle_mass]))
file.write('central mass:{}\n'.format([body_mass]))
file.write('forwarding rate:{}\n'.format(forwarding_rate))
file.write('one action time:{}\n'.format(one_action_time))
file.write('sleep time:{}\n'.format(sleep_time))
file.write('\n')

def force_X(theta, time, length, mass):
    theta = math.radians(theta)
    force = ((length*mass)*(theta/time)/2)*math.cos(theta)
    return force

def force_Y2(theta, time, length, mass):
    unit = 1
    if not theta==0:
        unit = theta/abs(theta)

    theta = math.radians(abs(theta))
    force = (((length*mass)*(theta/time)/2)*math.sin(theta))*unit*(-1)
    return force

def force_Y(theta, time, length, mass):
    unit = 1
    if not theta==0:
        unit = theta/abs(theta)

    theta = math.radians(theta)
    force = (((length*mass)*(theta/time)/2)*math.sin(theta))*unit
    return force

def just_force_Y(theta, time, length, mass):
    theta = math.radians(theta)
    force = -(((length * mass) * (theta / time) / 2) * math.sin(theta))
    return force

def deg_to_X23(x_start_point, r, deg):
    deg = math.radians(deg)
    X = -(math.sin(deg)*r)+x_start_point
    return X

def deg_to_X14(x_start_point, r, deg):
    deg = math.radians(deg)
    X = (math.sin(deg)*r)+x_start_point
    return X

def except_deg_to_X23(x_start_point, r, deg):
    unit = 1
    if deg != 0:
        unit = int(deg/abs(deg))
    deg = math.radians(deg)
    X = -(math.sin(deg)*r*unit)+x_start_point
    return X

def circle_eq(x_start_point, y_start_point, x, r):
    y = (((r**2)-((x-x_start_point)**2))**0.5)+y_start_point
    return y


def leg_center(leg_num, ancle, knee, upper, ancle_mass, knee_mass, upper_mass, degree_knee, degree_upper, knee_time, upper_time):
    deg_change = [0, 0]

    ancle_y = ancle[1]
    ancle_z = ancle[2]
    knee_y = knee[1]
    knee_z = knee[2]
    upper_y = upper[1]
    upper_z = upper[2]

    '''plt.plot([ancle_z], [ancle_y], 'go')
    plt.plot([knee_z], [knee_y], 'bo')
    plt.plot([upper_z], [upper_y], 'ro')
    plt.plot([upper_z, knee_z], [upper_y, knee_y], 'r')
    plt.plot([knee_z, ancle_z], [knee_y, ancle_y], 'r')'''

    standard_weight_center_Z = ((ancle_z*ancle_mass*g)+(knee_z*knee_mass*g)+(upper_z*upper_mass*g))/(ancle_mass*g+knee_mass*g+upper_mass*g)
    standard_weight_center_Y = ((ancle_y*ancle_mass*g)+(knee_y*knee_mass*g)+(upper_y*upper_mass*g))/(ancle_mass*g+knee_mass*g+upper_mass*g)
    standard_weight_center = (standard_weight_center_Z, standard_weight_center_Y)
    #print('standard:', standard_weight_center)

    #plt.plot([standard_weight_center_Z], [standard_weight_center_Y], 'ro')

    small_weight_center_z = ((ancle_z*ancle_mass*g)+(knee_z*knee_mass*g))/(ancle_mass*g+knee_mass*g)
    small_weight_center_y = ((ancle_y*ancle_mass*g)+(knee_y*knee_mass*g))/(ancle_mass*g+knee_mass*g)
    small_weight_center = (small_weight_center_z, small_weight_center_y)
    #print(small_weight_center)

    whole_time = max([knee_time, upper_time])
    whole_time_unit = whole_time/sleep_time

    knee_time_unit = knee_time/sleep_time
    upper_time_unit = upper_time/sleep_time

    devided_knee_degree_unit = degree_knee/knee_time_unit
    devided_knee_degree = []
    for i in range(int(knee_time_unit)):
        devided_knee_degree.append(devided_knee_degree_unit)
    for i in range(int(whole_time_unit-knee_time_unit)):
        devided_knee_degree.append(0.0)

    devided_upper_degree_unit = degree_upper/upper_time_unit
    devided_upper_degree = []
    for i in range(int(upper_time_unit)):
        devided_upper_degree.append(devided_upper_degree_unit)
    for i in range(int(whole_time_unit-upper_time_unit)):
        devided_upper_degree.append(0.0)

    #print(devided_knee_degree)
    #print(devided_upper_degree)
    #print(len(devided_upper_degree), len(devided_knee_degree))

    a, b = 0, 0
    if leg_num == 'A':
        a, b = 2, 3
    elif leg_num == 'B':
        a, b = 8, 9

    first_knee_Z = upper_z-knee_z
    first_knee_Y = upper_y-knee_y
    state = 1
    if first_knee_Z<0:
        state = -1
    tri = first_knee_Z/abs(first_knee_Y)
    match = 0
    count = 0
    original_deg_upper = 0
    while match==0:
        degtouse = count*0.01
        if tri-math.tan(math.radians(degtouse))<0.01:
            match = 1
            original_deg_upper = degtouse
        count = count+1
    original_deg_upper = original_deg_upper*state
    #print('knee original:', original_deg_upper)

    first_ancle_Z = knee_z-ancle_z
    first_ancle_Y = knee_y-ancle_y
    tri = abs(first_ancle_Z)/abs(first_ancle_Y)
    match = 0
    count = 0
    original_deg_knee = 0
    while match == 0:
        degtouse = count*0.01
        if tri-math.tan(math.radians(degtouse))<0.01:
            match = 1
            original_deg_knee = degtouse
        count = count+1
    #print('ancle original:', original_deg_knee)

    whole_weight_center_movements = []
    for i in range(int(whole_time_unit)):
        #global knee_z, knee_y, ancle_z, ancle_y
        deg_change[0] = deg_change[0]+devided_knee_degree[i]
        deg_change[1] = deg_change[1]+devided_upper_degree[i]
        distance_upper_to_weight_center = distance_2D((upper_z, upper_y), small_weight_center)
        #print(distance_upper_to_weight_center)

        weight_center_Z = ((ancle_z*(ancle_mass*g-force_X(deg_change[0], (i+1)*0.01, below_leg, ancle_mass)))+(knee_z*(knee_mass*g-force_X(deg_change[1], (i+1)*0.01, distance_upper_to_weight_center, ancle_mass+knee_mass)))+(upper_z*upper_mass*g))/(ancle_mass*g+knee_mass*g+upper_mass*g)
        #final formula!!!! check out forceX and forceY functions.


        #print('leg:', upper_leg, distance_upper_to_weight_center)
        #print('force result:', force_X(deg_change[0], (i+1)*0.01, below_leg, ancle_mass), force_X(deg_change[1], (i+1)*0.01, distance_upper_to_weight_center, ancle_mass+knee_mass))
        #print('ancle:', '({}, {})'.format(ancle_z, ancle_y))
        #print('knee:({}, {})'.format(knee_z, knee_y))
        #print('upper:({}, {})'.format(upper_z, upper_y))

        weight_center_Y = (ancle_y*(ancle_mass*g-force_Y(deg_change[0], (i+1)*0.01, below_leg, ancle_mass))+knee_y*(knee_mass*g-force_Y(deg_change[1], (i+1)*0.01, distance_upper_to_weight_center, ancle_mass+knee_mass))+upper_y*upper_mass*g)/(ancle_mass*g+knee_mass*g+upper_mass*g)

        #weight_center_Z = round(weight_center_Z, 3)
        #weight_center_Y = round(weight_center_Y, 3)
        weight_center = (weight_center_Z, weight_center_Y)
        #print('weight_center:', weight_center)


        old_deg[a] = old_deg[a]+deg_change[0]
        old_deg[b] = old_deg[b]+deg_change[1]

        #knee_z1 = deg_to_X23(upper_z, upper_leg, original_deg_upper-deg_change[1])
        #knee_y1 = (((upper_leg**2)-(knee_z1-upper_z)**2)**0.5)*(-1)+upper_y

        #ancle_z1 = deg_to_X14(knee_z1, below_leg, original_deg_knee+deg_change[0])
        #ancle_y1 = (((below_leg**2)-(ancle_z1-knee_z1)**2)**0.5)*(-1)+knee_y1

        knee_z = deg_to_X23(upper_z, upper_leg, original_deg_upper + deg_change[1])
        knee_y = (((upper_leg ** 2) - (knee_z - upper_z) ** 2) ** 0.5) * (-1) + upper_y

        ancle_z = deg_to_X14(knee_z, below_leg, original_deg_knee + deg_change[0]-deg_change[1])
        ancle_y = (((below_leg ** 2) - (ancle_z - knee_z) ** 2) ** 0.5) * (-1) + knee_y

        small_weight_center_z = ((ancle_z * (ancle_mass*g-force_X(deg_change[0], (i+1)*0.01, below_leg, ancle_mass))) + (knee_z * (knee_mass*g-force_X(deg_change[1], (i+1)*0.01, distance_upper_to_weight_center, ancle_mass+knee_mass)))) / (ancle_mass*g + knee_mass*g)
        small_weight_center_y = ((ancle_y * (ancle_mass*g-force_Y(deg_change[0], (i+1)*0.01, below_leg, ancle_mass))) + (knee_y * (knee_mass*g-force_Y(deg_change[1], (i+1)*0.01, distance_upper_to_weight_center, ancle_mass+knee_mass)))) / (ancle_mass*g + knee_mass*g)
        small_weight_center = (small_weight_center_z, small_weight_center_y)
        #print(small_weight_center)

        '''print(deg_change)
        print('ancle:', ancle_z, ancle_y)
        print('knee:', knee_z, knee_y)
        print(upper_z, upper_y)
        print('======')
        plt.plot([ancle_z], [ancle_y], 'go')
        plt.plot([knee_z], [knee_y], 'bo')
        plt.plot([upper_z], [upper_y], 'ro')
        plt.plot([upper_z, knee_z], [upper_y, knee_y], 'r')
        plt.plot([knee_z, ancle_z], [knee_y, ancle_y], 'r')
        plt.plot([weight_center_Z], [weight_center_Y], 'bo')'''

        #print('distance:', distance_2D((upper_z, upper_y), (knee_z, knee_y)))
        #print('distance:', distance_2D((knee_z, knee_y), (ancle_z, ancle_y)))

        just_weight_center_Z = ((ancle_z * ancle_mass * g) + (knee_z * knee_mass * g) + (upper_z * upper_mass * g)) / (ancle_mass * g + knee_mass * g + upper_mass * g)
        just_weight_center_Y = ((ancle_y * ancle_mass * g) + (knee_y * knee_mass * g) + (upper_y * upper_mass * g)) / (ancle_mass * g + knee_mass * g + upper_mass * g)
        #plt.plot([just_weight_center_Z], [just_weight_center_Y], 'go')


        whole_weight_center_movements.append(weight_center)

    #plt.show()
    #print(whole_weight_center_movements)

    return whole_weight_center_movements, (ancle[0], ancle_y, ancle_z), (knee[0], knee_y, knee_z)

def reverse_leg_center(leg_num, ancle, knee, upper, ancle_mass, knee_mass, upper_mass, degree_knee, degree_upper, knee_time, upper_time):
    deg_change = [0, 0]

    ancle_y = ancle[1]
    ancle_z = ancle[2]
    knee_y = knee[1]
    knee_z = knee[2]
    upper_y = upper[1]
    upper_z = upper[2]

    '''plt.plot([ancle_z], [ancle_y], 'go')
    plt.plot([knee_z], [knee_y], 'bo')
    plt.plot([upper_z], [upper_y], 'ro')
    plt.plot([upper_z, knee_z], [upper_y, knee_y], 'r')
    plt.plot([knee_z, ancle_z], [knee_y, ancle_y], 'r')'''

    standard_weight_center_Z = ((ancle_z*ancle_mass*g)+(knee_z*knee_mass*g)+(upper_z*upper_mass*g))/(ancle_mass*g+knee_mass*g+upper_mass*g)
    standard_weight_center_Y = ((ancle_y*ancle_mass*g)+(knee_y*knee_mass*g)+(upper_y*upper_mass*g))/(ancle_mass*g+knee_mass*g+upper_mass*g)
    standard_weight_center = (standard_weight_center_Z, standard_weight_center_Y)
    #print('standard:', standard_weight_center)

    #plt.plot([standard_weight_center_Z], [standard_weight_center_Y], 'ro')

    small_weight_center_z = ((upper_z*upper_mass*g)+(knee_z*knee_mass*g))/(upper_mass*g+knee_mass*g)
    small_weight_center_y = ((upper_y*upper_mass*g)+(knee_y*knee_mass*g))/(upper_mass*g+knee_mass*g)
    small_weight_center = (small_weight_center_z, small_weight_center_y)
    #print(small_weight_center)

    whole_time = max([knee_time, upper_time])
    whole_time_unit = whole_time/sleep_time

    knee_time_unit = knee_time/sleep_time
    upper_time_unit = upper_time/sleep_time

    devided_knee_degree_unit = degree_knee/knee_time_unit
    devided_knee_degree = []
    for i in range(int(knee_time_unit)):
        devided_knee_degree.append(devided_knee_degree_unit)
    for i in range(int(whole_time_unit-knee_time_unit)):
        devided_knee_degree.append(0.0)

    devided_upper_degree_unit = degree_upper/upper_time_unit
    devided_upper_degree = []
    for i in range(int(upper_time_unit)):
        devided_upper_degree.append(devided_upper_degree_unit)
    for i in range(int(whole_time_unit-upper_time_unit)):
        devided_upper_degree.append(0.0)

    #print(devided_knee_degree)
    #print(devided_upper_degree)
    #print(len(devided_upper_degree), len(devided_knee_degree))

    a, b = 0, 0
    if leg_num == 'A':
        a, b = 2, 3
    elif leg_num == 'B':
        a, b = 8, 9

    first_knee_Z = knee_z-ancle_z
    first_knee_Y = knee_y-ancle_y
    state = 1
    if first_knee_Z<0:
        state = -1
    tri = abs(first_knee_Z)/abs(first_knee_Y)
    match = 0
    count = 0
    original_deg_ancle = 0
    while match==0:
        degtouse = count*0.01
        if tri-math.tan(math.radians(degtouse))<0.01:
            match = 1
            original_deg_ancle = degtouse
        count = count+1
    original_deg_ancle = original_deg_ancle*state*(-1)
    #print('knee original:', original_deg_ancle)

    first_upper_Z = upper_z-knee_z
    first_upper_Y = upper_y-knee_y
    tri = abs(first_upper_Z)/abs(first_upper_Y)
    match = 0
    count = 0
    original_deg_knee = 0
    while match == 0:
        degtouse = count*0.01
        if tri-math.tan(math.radians(degtouse))<0.01:
            match = 1
            original_deg_knee = degtouse
        count = count+1
    #print('ancle original:', original_deg_knee)

    whole_weight_center_movements = []
    for i in range(int(whole_time_unit)):
        #global knee_z, knee_y, ancle_z, ancle_y
        deg_change[0] = deg_change[0]+devided_knee_degree[i]
        deg_change[1] = deg_change[1]+devided_upper_degree[i]
        distance_ancle_to_weight_center = distance_2D((ancle_z, ancle_y), small_weight_center)
        #print(distance_upper_to_weight_center)

        weight_center_Z = ((ancle_z*ancle_mass*g)+(knee_z*(knee_mass*g-force_X(deg_change[0], (i+1)*0.01, distance_ancle_to_weight_center, ancle_mass+knee_mass)))+(upper_z*(upper_mass*g-force_X(deg_change[1], (i+1)*0.01, upper_leg, upper_mass))))/(ancle_mass*g+knee_mass*g+upper_mass*g)
        #final formula!!!! check out forceX and forceY functions.


        #print('leg:', upper_leg, distance_upper_to_weight_center)
        #print('force result:', force_X(deg_change[0], (i+1)*0.01, below_leg, ancle_mass), force_X(deg_change[1], (i+1)*0.01, distance_upper_to_weight_center, ancle_mass+knee_mass))
        #print('ancle:', '({}, {})'.format(ancle_z, ancle_y))
        #print('knee:({}, {})'.format(knee_z, knee_y))
        #print('upper:({}, {})'.format(upper_z, upper_y))

        weight_center_Y = ((ancle_y*ancle_mass*g)+(knee_y*(knee_mass*g-force_Y2(deg_change[0], (i+1)*0.01, distance_ancle_to_weight_center, ancle_mass+knee_mass)))+(upper_y*(upper_mass*g-force_Y2(deg_change[1], (i+1)*0.01, upper_leg, upper_mass))))/(ancle_mass*g+knee_mass*g+upper_mass*g)

        #weight_center_Z = round(weight_center_Z, 3)
        #weight_center_Y = round(weight_center_Y, 3)
        weight_center = (weight_center_Z, weight_center_Y)
        #print('weight_center:', weight_center)


        old_deg[a] = old_deg[a]+deg_change[0]
        old_deg[b] = old_deg[b]+deg_change[1]

        #knee_z1 = deg_to_X23(upper_z, upper_leg, original_deg_upper-deg_change[1])
        #knee_y1 = (((upper_leg**2)-(knee_z1-upper_z)**2)**0.5)*(-1)+upper_y

        #ancle_z1 = deg_to_X14(knee_z1, below_leg, original_deg_knee+deg_change[0])
        #ancle_y1 = (((below_leg**2)-(ancle_z1-knee_z1)**2)**0.5)*(-1)+knee_y1

        knee_z = deg_to_X23(ancle_z, below_leg, original_deg_ancle + deg_change[0])
        knee_y = (((below_leg ** 2) - (knee_z - ancle_z) ** 2) ** 0.5) + ancle_y

        upper_z = deg_to_X14(knee_z, upper_leg, original_deg_knee - deg_change[0]+deg_change[1])
        upper_y = (((upper_leg ** 2) - (upper_z - knee_z) ** 2) ** 0.5) + knee_y

        small_weight_center_z = ((upper_z * (upper_mass*g-force_X(deg_change[1], (i+1)*0.01, upper_leg, upper_mass))) + (knee_z * (knee_mass*g-force_X(deg_change[0], (i+1)*0.01, distance_ancle_to_weight_center, upper_mass+knee_mass)))) / (upper_mass*g + knee_mass*g)
        small_weight_center_y = ((upper_y * (upper_mass*g-force_Y2(deg_change[1], (i+1)*0.01, upper_leg, upper_mass))) + (knee_y * (knee_mass*g-force_Y2(deg_change[0], (i+1)*0.01, distance_ancle_to_weight_center, upper_mass+knee_mass)))) / (upper_mass*g + knee_mass*g)
        small_weight_center = (small_weight_center_z, small_weight_center_y)
        #print(small_weight_center)

        '''print(deg_change)
        print('ancle:', ancle_z, ancle_y)
        print('knee:', knee_z, knee_y)
        print(upper_z, upper_y)
        print('======')
        plt.plot([ancle_z], [ancle_y], 'go')
        plt.plot([knee_z], [knee_y], 'bo')
        plt.plot([upper_z], [upper_y], 'ro')
        plt.plot([upper_z, knee_z], [upper_y, knee_y], 'r')
        plt.plot([knee_z, ancle_z], [knee_y, ancle_y], 'r')
        plt.plot([weight_center_Z], [weight_center_Y], 'bo')'''

        #print('distance:', distance_2D((upper_z, upper_y), (knee_z, knee_y)))
        #print('distance:', distance_2D((knee_z, knee_y), (ancle_z, ancle_y)))

        just_weight_center_Z = ((ancle_z * ancle_mass * g) + (knee_z * knee_mass * g) + (upper_z * upper_mass * g)) / (ancle_mass * g + knee_mass * g + upper_mass * g)
        just_weight_center_Y = ((ancle_y * ancle_mass * g) + (knee_y * knee_mass * g) + (upper_y * upper_mass * g)) / (ancle_mass * g + knee_mass * g + upper_mass * g)
        #plt.plot([just_weight_center_Z], [just_weight_center_Y], 'go')


        whole_weight_center_movements.append(weight_center)

    #plt.show()
    #print(whole_weight_center_movements)

    return whole_weight_center_movements, (upper[0], upper_y, upper_z), (knee[0], knee_y, knee_z)

def weight_center_3D(list_of_particles, masses):
    xs = []
    ys = []
    zs = []
    for particle in list_of_particles:
        xs.append(particle[0])
        ys.append(particle[1])
        zs.append(particle[2])

    x_sum, y_sum, z_sum = 0, 0, 0
    for i in range(len(xs)):
        x_sum = x_sum + xs[i]*masses[i]
        y_sum = y_sum + ys[i]*masses[i]
        z_sum = z_sum + zs[i]*masses[i]

    mass_sum = 0
    for i in masses:
        mass_sum = mass_sum + i

    x_sum = x_sum/mass_sum
    y_sum = y_sum/mass_sum
    z_sum = z_sum/mass_sum
    final_list = (x_sum, y_sum, z_sum)
    return final_list

def side_moves(right_leg, left_leg, right_foot, left_foot, right_mass, left_mass, central_mass, central_particle, ancle_theta, first_ancle_theta, time, wait_time):
    global right_x, right_y, left_x, left_y, central_x, central_y
    right_x, right_y = right_leg[0], right_leg[1]
    right_foot_x, right_foot_y = right_foot[0], right_foot[1]
    left_x, left_y = left_leg[0], left_leg[1]
    left_foot_x, left_foot_y = left_foot[0], left_foot[1]
    central_x, central_y = central_particle[0], central_particle[1]
    #const_right_x, const_right_y, const_left_x, const_left_y, const_central_x, const_central_y = right_x, right_y, left_x, left_y, central_x, central_y
    const_central_x, const_central_y = (right_foot_x+left_foot_x)/2, (right_foot_y+left_foot_y)/2
    central_to_foot_middle_distance = distance_2D((const_central_x, const_central_y), (central_x, central_y))

    right_leg_len = distance_2D((right_x, right_y), (right_foot_x, right_foot_y))
    left_leg_len = distance_2D((left_x, left_y), (left_foot_x, left_foot_y))

    '''plt.plot([right_x], [right_y], 'go')
    plt.plot([left_x], [left_y], 'go')
    plt.plot([central_x], [central_y], 'go')'''

    steps = int(time/0.01)
    unit = ancle_theta/steps
    deg_change_list = []
    for i in range(steps):
        deg_change_list.append(unit)
    for i in range(int(wait_time/0.01)):
        deg_change_list.append(0)

    change_sum = first_ancle_theta
    pure_change = 0
    whole_weight_center_changes = []
    for i, change in enumerate(deg_change_list):
        change_sum = change_sum + change
        pure_change = pure_change + change
        if ancle_theta<0:
            right_x = except_deg_to_X23(x_start_point=right_foot_x, r=right_leg_len, deg=change_sum)
            right_y = circle_eq(x_start_point=right_foot_x, y_start_point=right_foot_y, x=right_x, r=right_leg_len)
            left_x = except_deg_to_X23(x_start_point=left_foot_x, r=left_leg_len, deg=change_sum)
            left_y = circle_eq(x_start_point=left_foot_x, y_start_point=left_foot_y, x=left_x, r=left_leg_len)
            central_x = except_deg_to_X23(x_start_point=const_central_x, r=central_to_foot_middle_distance, deg=change_sum)
            central_y = circle_eq(x_start_point=const_central_x, y_start_point=const_central_y, r=central_to_foot_middle_distance, x=central_x)

            moved_weight_center_X = ((right_x*((right_mass*g)-force_X(pure_change, (i+1)*0.01, right_leg_len, right_mass)))+(left_x*((left_mass*g)-force_X(pure_change, (i+1)*0.01, left_leg_len, left_mass)))+(central_x*((central_mass*g)-force_X(pure_change, (i+1)*0.01, central_to_foot_middle_distance, central_mass))))/((right_mass+left_mass+central_mass)*g)
            moved_weight_center_Y = ((right_y*((right_mass*g)-just_force_Y(pure_change, (i+1)*0.01, right_leg_len, right_mass)))+(left_y*((left_mass*g)-just_force_Y(pure_change, (i+1)*0.01, left_leg_len, left_mass)))+(central_y*((central_mass*g)-just_force_Y(pure_change, (i+1)*0.01, central_to_foot_middle_distance, central_mass))))/((right_mass+left_mass+central_mass)*g)
            whole_weight_center_changes.append([moved_weight_center_X, moved_weight_center_Y])
            '''plt.plot([moved_weight_center_X], [moved_weight_center_Y], 'go')

            plt.plot([right_x], [right_y], 'ro')
            plt.plot([left_x], [left_y], 'ro')
            plt.plot([central_x], [central_y], 'bo')'''
        else:
            right_x = deg_to_X14(x_start_point=right_foot_x, r=right_leg_len, deg=change_sum)
            right_y = circle_eq(x_start_point=right_foot_x, y_start_point=right_foot_y, r=right_leg_len, x=right_x)
            left_x = deg_to_X14(x_start_point=left_foot_x, r=left_leg_len, deg=change_sum)
            left_y = circle_eq(x_start_point=left_foot_x, y_start_point=left_foot_y, r=left_leg_len, x=left_x)
            central_x = deg_to_X14(x_start_point=const_central_x, r=central_to_foot_middle_distance, deg=change_sum)
            central_y = circle_eq(x_start_point=const_central_x, y_start_point=const_central_y, r=central_to_foot_middle_distance, x=central_x)

            moved_weight_center_X = ((right_x * ((right_mass * g) - force_X(pure_change, (i + 1) * 0.01, right_leg_len, right_mass))) + (
                                                 left_x * ((left_mass * g) - force_X(pure_change, (i + 1) * 0.01,
                                                                                     left_leg_len, left_mass))) + (
                                                 central_x * ((central_mass * g) - force_X(pure_change, (i + 1) * 0.01,
                                                                                           central_to_foot_middle_distance,
                                                                                           central_mass)))) / (
                                                (right_mass + left_mass + central_mass) * g)
            moved_weight_center_Y = ((right_y * (
                        (right_mass * g) - force_Y(pure_change, (i + 1) * 0.01, right_leg_len, right_mass))) + (
                                                 left_y * ((left_mass * g) - just_force_Y(pure_change, (i + 1) * 0.01,
                                                                                     left_leg_len, left_mass))) + (
                                                 central_y * ((central_mass * g) - just_force_Y(pure_change, (i + 1) * 0.01,
                                                                                           central_to_foot_middle_distance,
                                                                                           central_mass)))) / (
                                                (right_mass + left_mass + central_mass) * g)
            #plt.plot([moved_weight_center_X], [moved_weight_center_Y], 'go')
            whole_weight_center_changes.append([moved_weight_center_X, moved_weight_center_Y])

            '''plt.plot([right_x], [right_y], 'ro')
            plt.plot([left_x], [left_y], 'ro')
            plt.plot([central_x], [central_y], 'bo')'''

    #plt.show()

    return whole_weight_center_changes, (right_x, right_y, right_leg[2]), (left_x, left_y, left_leg[2]), (central_x, central_y, central_particle[2])

def equation2D(factors):
    [a, b, c] = factors
    r1 = (-b+(b**2-4*(a*c))**0.5)/(2*a)
    r2 = (-b - (b ** 2 - 4 * (a * c)) ** 0.5) / (2 * a)
    return [r1, r2]

def line_generator(a, b):
    ax, ay = a[0], a[1]
    bx, by = b[0], b[1]
    how_steep = (by-ay)/(bx-ax)
    rest = ay-how_steep*ax
    def line(x):
        y = how_steep*x+rest
        return y

    return line, how_steep, rest

def sin_reverser(result):
    deg = 0
    for i in range(90*1000000):
        middle = math.radians(deg)
        compare = math.sin(middle)
        if abs(compare-result)<0.0001:
            return deg
        deg = deg + i/1000000

def degree_finder(A, B, C):
    Az, Ay = A[2], A[1]
    Bz, By = B[2], B[1]
    Cz, Cy = C[2], C[1]

    y_differential_A2B = abs(Ay-By)
    y_differential_B2C = abs(By-Cy)
    small_A = (Az, Ay)
    small_B = (Bz, By)
    small_C = (Cz, Cy)
    result1 = (y_differential_A2B)/distance_2D(small_A, small_B)
    result2 = (y_differential_B2C)/distance_2D(small_B, small_C)
    theta1 = sin_reverser(result1)
    theta2 = sin_reverser(result2)
    theta = theta1+theta2

    if Bz<Az and Bz<Cz:
        side = 'left'
    else:
        side = 'right'

    return theta, side

def get_leg_steep(top, bottom):
    topx, topy = top[0], top[1]
    bottomx, bottomy = bottom[0], bottom[1]
    tuplizex = (topx, topy)
    tuplizey = (bottomx, bottomy)
    length = distance_2D(tuplizey, tuplizex)
    width_differential = abs(topx-bottomx)
    sin_theta = width_differential/length
    theta = sin_reverser(sin_theta)

    if topx>bottomx:
        side = 'right'
    else:
        side = 'left'

    return theta, side

def get_knee_points(body_fix, foot, lower_leg, upper_leg):
    foot_y, foot_z = foot[1], foot[2]
    body_y, body_z = distance_2D((foot[0], foot[1]), (body_fix[0], body_fix[1]))+foot_y, body_fix[2]

    z_range_start = foot_z-lower_leg
    z_range_end = foot_z+lower_leg

    repeat_count = round(z_range_end-z_range_start, 2)

    knee_positions = []

    for i in range(int(repeat_count*10000)):
        x = z_range_start+i/10000
        eq = (lower_leg**2-(x-foot_z)**2)**0.5+foot_y+(upper_leg**2-(x-body_z)**2)**0.5-body_y

        if type(eq)==float:
            if round(eq, 3)==0:
                y = (lower_leg**2-(x-foot_z)**2)**0.5+foot_y
                knee_positions.append((y, x))

    '''for i in knee_positions:
        print('body', distance_2D(i, (body_y, body_z)))
        print('foot', distance_2D(i, (foot_y, foot_z)))'''


    knee_positions_z = []
    knee_positions_y = []
    for i in knee_positions:
        knee_positions_z.append(i[1])
        knee_positions_y.append(i[0])

    selected = []
    if len(knee_positions)>2:
        selectedA = max(knee_positions_z)
        selectedB = min(knee_positions_z)
        selectedA_num = knee_positions_z.index(selectedA)
        selectedB_num = knee_positions_z.index(selectedB)
        selectedYA = knee_positions_y[selectedA_num]
        selectedYB = knee_positions_y[selectedB_num]
        selected.append((selectedYB, selectedB))
        selected.append((selectedYA, selectedA))

    if foot[0]!=body_fix[0]:
        add_x = []
        turned_degree, side = get_leg_steep(body_fix, foot)
        if len(knee_positions)>2:
            if side == 'right':
                for i in selected:
                    y = i[0]
                    x = deg_to_X14(foot[0], y, turned_degree)
                    adjust_y = (y**2-(x-foot[0])**2)**0.5+foot[1]
                    new = (x, adjust_y, i[1])
                    add_x.append(new)
            else:
                for i in selected:
                    y = i[0]
                    x = deg_to_X23(foot[0], y, turned_degree)
                    adjust_y = (y**2-(x-foot[0])**2)**0.5+foot[1]
                    new = (x, adjust_y, i[1])
                    add_x.append(new)
        else:
            if side == 'right':
                for i in knee_positions:
                    y = i[0]
                    x = deg_to_X14(foot[0], y, turned_degree)
                    adjust_y = (y ** 2 - (x - foot[0]) ** 2) ** 0.5 + foot[1]
                    new = (x, adjust_y, i[1])
                    add_x.append(new)
            else:
                for i in knee_positions:
                    y = i[0]
                    x = deg_to_X23(foot[0], y, turned_degree)
                    adjust_y = (y ** 2 - (x - foot[0]) ** 2) ** 0.5 + foot[1]
                    new = (x, adjust_y, i[1])
                    add_x.append(new)

        return add_x

    else:
        add_x = []
        if len(knee_positions)>2:
            for i in selected:
                new = (foot[0], i[0], i[1])
                add_x.append(new)
        else:
            for i in knee_positions:
                new = (foot[0], i[0], i[1])
                add_x.append(new)

        return add_x

def leg_step_generator(how_far, time, body_fix, foot_origin, lower_leg, upper_leg):
    foot_movement = []
    how_many_steps = int(time/0.01)
    one_step = how_far/how_many_steps

    foot_z, foot_y = foot_origin[2], foot_origin[1]

    for i in range(how_many_steps+1):
        adjust_foot_y = ((how_far**2)/4-((foot_z+one_step*i)-(foot_z+how_far/2))**2)**0.5+foot_y
        tuplize = (adjust_foot_y, foot_z+one_step*i)
        foot_movement.append(tuplize)

    x_find = []
    if body_fix[0]!=foot_origin[0]:
        turned_degree, side = get_leg_steep(body_fix, foot_origin)
        if side == 'right':
            for i in foot_movement:
                y = i[0]-foot_origin[1]
                x = deg_to_X14(x_start_point=foot_origin[0], r=y, deg=turned_degree)
                adjust_y = (y**2-(x-foot_origin[0])**2)**0.5+foot_origin[1]
                new = (x, adjust_y, i[1])
                x_find.append(new)
        else:
            for i in foot_movement:
                y = i[0]
                x = deg_to_X23(foot_origin[0], y-foot_origin[1], turned_degree)
                adjust_y = ((y-foot_origin[1])**2-(x-foot_origin[0])**2)**0.5+foot_origin[1]
                new = (x, adjust_y, i[1])
                x_find.append(new)

    else:
        for i in foot_movement:
            new = (body_fix[0], i[0], i[1])
            x_find.append(new)

    knee_points = []
    for i in x_find:
        knee_point = get_knee_points(body_fix, i, lower_leg, upper_leg)
        for a in knee_point:
            degree, side = degree_finder(body_fix, a, i)
            if side == 'left':
                knee_points.append(a)

    final_result = []
    for i in knee_points:
        new = (body_fix, i, x_find[knee_points.index(i)])
        final_result.append(new)

    degs = []
    for i in final_result:
        target = (foot_origin[0], foot_origin[1])
        foot_tuple = (i[2][0], i[2][1])
        knee_tuple = (i[1][0], i[1][1])
        body_tuple = (i[0][0], i[0][1])

        new_foot_y = distance_2D(target, foot_tuple)+foot_origin[1]
        new_knee_y = distance_2D(target, knee_tuple)+foot_origin[1]
        new_body_y = distance_2D(target, body_tuple)+foot_origin[1]

        new_foot = (foot_origin[0], new_foot_y, i[2][2])
        new_knee = (foot_origin[0], new_knee_y, i[1][2])
        new_body = (foot_origin[0], new_body_y, i[0][2])

        for_steep_body = (new_body[2], new_body[1])
        for_steep_knee = (new_knee[2], new_knee[1])
        for_steep_foot = (new_foot[2], new_foot[1])

        upper_steep, side1 = get_leg_steep(for_steep_body, for_steep_knee)
        lower_steep, side2 = get_leg_steep(for_steep_knee, for_steep_foot)

        degs.append((upper_steep, lower_steep))

    return final_result, degs


#leg_center(leg_num='A', ancle=ancle_A_first_weight_point, knee=knee_A_weight_point, upper=back_forth_A_weight_point, ancle_mass=below_leg_mass, knee_mass=upper_leg_mass, upper_mass=body_mass, degree_upper=-30, degree_knee=0, upper_time=0.01, knee_time=0.01)
#reverse_leg_center(leg_num='A', ancle=ancle_A_first_weight_point, knee=knee_A_weight_point, upper=back_forth_A_weight_point, ancle_mass=below_leg_mass, knee_mass=upper_leg_mass, upper_mass=body_mass, degree_upper=0, degree_knee=-10, upper_time=0.01, knee_time=0.01)


#print(weight_center_3D([(1, 2, 3), (1, 2, 3), (1, 2, 3)], [1, 1, 1]))

#print(side_moves((1, 2, 1), (4, 2, 1), (1, 0, 2), (4, 0, 2), 3, 3, 3, (2.5, 5, 2), 0, 0, 0.01, wait_time=1))

#print(get_knee_points((1, 8, 3), (1, 0, 3), 5, 5))

def deg_curve(particle, shaft, deg):
    particlex = particle[0]
    particley = particle[1]
    shaftx, shafty = shaft[0], shaft[1]
    r = distance_2D((particlex, particley), (shaftx, shafty))
    adjust_x = math.sin(math.radians(deg))*r+shaftx
    adjust_y = (r**2-(adjust_x-shaftx)**2)**0.5+shafty
    return (adjust_x, adjust_y, particle[2])

upper = 0
knee = 1
leg1_step, leg1_degs = leg_step_generator(how_far=forwarding_rate, time=one_action_time, body_fix=leg1_upper_particle, foot_origin=leg1_ancle_particle, lower_leg=below_leg, upper_leg=upper_leg)
status = round(12.5)
guage()

#print(leg1_degs)
first_upper_deg_in_steps, first_knee_deg_in_steps = leg1_degs[0][0], leg1_degs[0][1]
latest_deg_upper, latest_deg_knee = first_upper_deg_in_steps, first_knee_deg_in_steps
deg_changes_upper, deg_changes_knee = [0.0], [0.0]
for i in leg1_degs:
    change_upper = i[upper]-latest_deg_upper
    deg_changes_upper.append(change_upper)
    latest_deg_upper = i[upper]

    change_knee = i[knee]-latest_deg_knee
    deg_changes_knee.append(change_knee)
    latest_deg_knee = i[knee]

#print(deg_changes_upper)
#print(deg_changes_knee)

leg1_weight_center_movements = []
for i in range(len(leg1_step)):
    leg1_weight_center, anclezy, kneezy = leg_center(leg_num='A', ancle=leg1_step[i][2], knee=leg1_step[i][1], upper=leg1_step[i][0], ancle_mass=leg1_ancle_mass, knee_mass=leg1_knee_mass, upper_mass=leg1_upper_mass, degree_upper=deg_changes_upper[i], degree_knee=deg_changes_knee[i], upper_time=sleep_time, knee_time=sleep_time)
    manipulated = (leg1_step[i][0][0], leg1_weight_center[0][1], leg1_weight_center[0][0])
    leg1_weight_center_movements.append(manipulated)

foot_arrived_position = leg1_step[-1]

still_leg_upper = leg2_upper_particle
still_leg_knee = leg2_knee_particle
still_leg_ancle = leg2_ancle_particle
still_leg_weight_center = weight_center_3D([still_leg_upper, still_leg_knee, still_leg_ancle], [leg2_upper_mass, leg2_knee_mass, leg2_ancle_mass])

body_weight_center = body_weight_particle
body_weight_center_mass = body_mass
target_position = leg2_ancle_particle

degrees = []   #body move final result!
pre_match = 0
degree_for_finding = 0
while pre_match == 0:
    #print(degree_for_finding)
    weight_center_changes, right_leg_change, left_leg_change, central_change = side_moves(right_leg=leg1_weight_center_movements[0],
                                                                                          left_leg=still_leg_weight_center,
                                                                                          right_foot=leg1_step[0][2],
                                                                                          left_foot=still_leg_ancle,
                                                                                          right_mass=leg1_upper_mass+leg1_knee_mass+leg1_ancle_mass, left_mass=leg2_upper_mass+leg2_knee_mass+leg2_ancle_mass,
                                                                                          central_mass=body_weight_center_mass,
                                                                                          central_particle=body_weight_center,
                                                                                          ancle_theta=degree_for_finding * 0.1,
                                                                                          first_ancle_theta=0,
                                                                                          time=sleep_time,
                                                                                          wait_time=3)
    if abs(weight_center_changes[-1][0]-target_position[0])<gap_allow:
        degrees.append(degree_for_finding*0.1)
        pre_match = 1
        #print('found', degree_for_finding*0.1)

    degree_for_finding += 1

right_leg_movement = []
left_leg_movement = []
central_movement = []
bowls = []
for a, i in enumerate(leg1_weight_center_movements):
    degree = 0
    #print(degree)
    match = 0
    bowl = []
    count = 100


    #print(count)
    while match == 0 and len(bowl) != count:
        last_deg = 0
        try:
            if degrees[-1] != 0:
                last_deg = degrees[-1]
        except:
            st = 'no last deg'
            #print('no last deg')
        adjust_right_leg = deg_curve(i, leg1_step[a][2], last_deg)
        adjust_left_leg = deg_curve(still_leg_weight_center, still_leg_ancle, last_deg)
        adjust_central = deg_curve(body_weight_center, ((leg1_step[a][2][0]+still_leg_ancle[0])/2, leg1_step[a][2][1]), last_deg)
        weight_center_changes, right_leg_change, left_leg_change, central_change = side_moves(right_leg=adjust_right_leg, left_leg=adjust_left_leg, right_foot=leg1_step[a][2], left_foot=still_leg_ancle, right_mass=leg1_upper_mass+leg1_knee_mass+leg1_ancle_mass, left_mass=leg2_upper_mass+leg2_knee_mass+leg2_ancle_mass, central_mass=body_weight_center_mass, central_particle=adjust_central, ancle_theta=degree*0.1, first_ancle_theta=last_deg, time=sleep_time*(a+1), wait_time=0)
        #print(weight_center_3D([i, still_leg_weight_center, body_weight_center], [5, 5, 3])[2])
        #print(last_deg)
        degree += 1
        '''if len(bowl)==0:
            target = 0
        else:
            target = bowl[-1]'''
        #print(target-degree*0.1, ',', target)
        #if abs(target-degree*0.1)>1:
        if abs(weight_center_changes[-1][0] - target_position[0]) < gap_allow*4:
            bowl.append(degree * 0.1 + last_deg)
            right_leg_movement.append(right_leg_change)
            left_leg_movement.append(left_leg_movement)
            central_movement.append(central_movement)
            if len(bowl) == count:
                match = 1

            if degree * 0.1 > 90:
                break


    #print(bowl)
    to_compare = degrees[0]
    bowl2 = []
    for b in bowl:
        thing = abs(to_compare-b)
        bowl2.append(thing)
    minimum_one = min(bowl2)
    min_index = bowl2.index(minimum_one)
    the_degree = bowl[min_index]
    degrees.append(the_degree)
    bowls.append(bowl)
    if (a+1)==(len(leg1_weight_center_movements)-1)/2:
        break


for i in range(len(degrees)-1):
    degrees.append(degrees[len(degrees)-2-2*i])

status = 18
guage()
#print('count', len(degrees))
#print(bowls)
#print(degrees)
body_forward = forwarding_rate
latest_leg1_upper, latest_leg1_knee, latest_leg1_ancle = leg1_step[-1]

target_upper_point_body_move = (latest_leg1_upper[0], latest_leg1_upper[1], latest_leg1_upper[2]+body_forward)
target_knee_body_move = get_knee_points(body_fix=target_upper_point_body_move, foot=latest_leg1_ancle, lower_leg=below_leg, upper_leg=upper_leg)
knee_point = 0
for i in target_knee_body_move:
    theta, side = degree_finder(target_upper_point_body_move, i, latest_leg1_ancle)
    if side == 'left':
        knee_point = i

for_steep_upper1 = (target_upper_point_body_move[2], target_upper_point_body_move[1])
for_steep_knee1 = (knee_point[2], knee_point[1])
for_steep_ancle1 = (latest_leg1_ancle[2], latest_leg1_ancle[1])
upper_deg_leg1, upper_side = get_leg_steep(for_steep_upper1, for_steep_knee1)
knee_deg_leg1, knee_side = get_leg_steep(for_steep_knee1, for_steep_ancle1)

for_steep_upper1_no_move = (latest_leg1_upper[2], latest_leg1_upper[1])
for_steep_knee1_no_move = (latest_leg1_knee[2], latest_leg1_knee[1])
for_steep_ancle1_no_move = (latest_leg1_ancle[2], latest_leg1_ancle[1])
upper_deg_leg1_no_move, no_move_upper_side = get_leg_steep(for_steep_upper1_no_move, for_steep_knee1_no_move)
knee_deg_leg1_no_move, no_move_knee_side = get_leg_steep(for_steep_knee1_no_move, for_steep_ancle1_no_move)

leg1_body_move_upper = upper_deg_leg1-upper_deg_leg1_no_move
leg1_body_move_knee = knee_deg_leg1-knee_deg_leg1_no_move

leg1_body_move_leg_movements = [leg1_body_move_upper, leg1_body_move_knee]

latest_leg2_upper, latest_leg2_knee, latest_leg2_ancle = still_leg_upper, still_leg_knee, still_leg_ancle

target_upper_point_body_move_2 = (latest_leg2_upper[0], latest_leg2_upper[1], latest_leg2_upper[2]+body_forward)
target_knee_body_move_2 = get_knee_points(body_fix=target_upper_point_body_move_2, foot=latest_leg2_ancle, lower_leg=below_leg, upper_leg=upper_leg)
knee_point2 = 0
for i in target_knee_body_move_2:
    theta, side = degree_finder(target_upper_point_body_move_2, i, latest_leg2_ancle)
    if side == 'left':
        knee_point2 = i

for_steep_upper2 = (target_upper_point_body_move_2[2], target_upper_point_body_move_2[1])
for_steep_knee2 = (knee_point2[2], knee_point2[1])
for_steep_ancle2 = (latest_leg2_ancle[2], latest_leg2_ancle[1])
upper_deg_leg2, upper_deg_side = get_leg_steep(for_steep_upper2, for_steep_knee2)
knee_deg_leg2, knee_deg_side = get_leg_steep(for_steep_knee2, for_steep_ancle2)

for_steep_upper2_no_move = (latest_leg2_upper[2], latest_leg2_upper[1])
for_steep_knee2_no_move = (latest_leg2_knee[2], latest_leg2_knee[1])
for_steep_ancle2_no_move = (latest_leg2_ancle[2], latest_leg2_ancle[1])
upper_deg_leg2_no_move, upper_deg_side2 = get_leg_steep(for_steep_upper2_no_move, for_steep_knee2_no_move)
knee_deg_leg2_no_move, knee_deg_side2 = get_leg_steep(for_steep_knee2_no_move, for_steep_ancle2_no_move)

leg2_body_move_upper = upper_deg_leg2-upper_deg_leg2_no_move
leg2_body_move_knee = knee_deg_leg2-knee_deg_leg2_no_move

leg2_body_move_leg_movements = [leg2_body_move_upper, leg2_body_move_knee]

leg1_body_move_weight_centers, leg1_upper, leg1_knee = reverse_leg_center(leg_num='A', ancle=latest_leg1_ancle, knee=latest_leg1_knee, upper=latest_leg1_upper, ancle_mass=leg1_ancle_mass, knee_mass=leg1_knee_mass, upper_mass=leg1_upper_mass, degree_knee=leg1_body_move_leg_movements[1], degree_upper=leg1_body_move_leg_movements[0], knee_time=one_action_time, upper_time=one_action_time)
leg2_body_move_weight_centers, leg2_upper, leg2_knee = reverse_leg_center(leg_num='B', ancle=latest_leg2_ancle, knee=latest_leg2_knee, upper=latest_leg2_upper, ancle_mass=leg2_ancle_mass, knee_mass=leg2_knee_mass, upper_mass=leg2_upper_mass, degree_knee=leg2_body_move_leg_movements[1], degree_upper=leg2_body_move_leg_movements[0], knee_time=one_action_time, upper_time=one_action_time)

one_body_forward_unit = body_forward/(one_action_time*100)
body_forward_current = 0
body_move_leg1_steps = [leg1_step[-1]]
for i in range(50):
    ancle1 = leg1_step[-1][2]
    upper1 = (leg1_step[-1][0][0], leg1_step[-1][0][1], leg1_step[-1][0][2]+((i+1)*one_body_forward_unit))
    knee_points = get_knee_points(upper1, ancle1, upper_leg=upper_leg, lower_leg=below_leg)
    correct_point = 0
    for a in knee_points:
        degree, side = degree_finder(upper1, a, ancle1)
        if side == 'left':
            correct_point = a
    to_return = [upper1, correct_point, ancle1]
    body_move_leg1_steps.append(to_return)
    body_forward_current += one_body_forward_unit

body_forward_current2 = 0
body_move_leg2_steps = [(latest_leg2_upper, latest_leg2_knee, latest_leg2_ancle)]
for i in range(50):
    ancle1 = latest_leg2_ancle
    upper1 = (latest_leg2_upper[0], latest_leg2_upper[1], latest_leg2_upper[2]+((i+1)*one_body_forward_unit))
    knee_points = get_knee_points(upper1, ancle1, upper_leg=upper_leg, lower_leg=below_leg)
    correct_point = 0
    for a in knee_points:
        degree, side = degree_finder(upper1, a, ancle1)
        if side == 'left':
            correct_point = a
    to_return = [upper1, correct_point, ancle1]
    body_move_leg2_steps.append(to_return)
    body_forward_current2 += one_body_forward_unit

status = 50
guage()
#print('body move 1:', body_move_leg1_steps)
#print('body move 2:', body_move_leg2_steps)

#print(leg1_upper)
central_particle_weight_centers = []
for i in range(int(one_action_time*100)):
    moved_amount = (body_forward/(one_action_time*100))*(i+1)
    adjusting_force = (body_weight_center_mass*moved_amount)/(2*((i+1)/100)**2)
    tuplize = (body_weight_center[2], body_weight_center[1])
    adjust_z = (tuplize[0]*(body_weight_center_mass*g-adjusting_force))/(body_weight_center_mass*g)
    adjusted = (body_weight_center[0], tuplize[1], adjust_z)
    central_particle_weight_centers.append(adjusted)

#print(central_particle_weight_centers[-1], leg1_body_move_weight_centers[-1], leg2_body_move_weight_centers[-1])

stablize_degree = 0
match_for_stablize = 0
stablized_degree = 0
manipulated1 = (latest_leg1_ancle[0], leg1_body_move_weight_centers[-1][1], leg1_body_move_weight_centers[-1][0])
manipulated2 = (latest_leg2_ancle[0], leg2_body_move_weight_centers[-1][1], leg2_body_move_weight_centers[-1][0])
while match_for_stablize == 0:
    weight_center_changes, right_leg_change, left_leg_change, central_change = side_moves(right_leg=manipulated1, left_leg=manipulated2, right_foot=latest_leg1_ancle, left_foot=latest_leg2_ancle, right_mass=leg1_upper_mass+leg1_knee_mass+leg1_ancle_mass, left_mass=leg2_upper_mass+leg2_knee_mass+leg2_ancle_mass, central_mass=body_weight_center_mass, central_particle=central_particle_weight_centers[-1], ancle_theta=-(stablize_degree*0.1), first_ancle_theta=0, time=sleep_time, wait_time=5)
    if abs(latest_leg1_ancle[0]-weight_center_changes[-1][0]) < gap_allow:
        match_for_stablize = 1
        stablized_degree = -stablize_degree*0.1

    stablize_degree += 1
    #print(stablize_degree)

#print('stablized:', stablized_degree)

one_step_unit = ((abs(stablized_degree)+abs(degrees[-1])))/(len(leg1_body_move_weight_centers))
last_deg2 = degrees[-1]
adjusted_body_move_leg1_steps = []
adjusted_body_move_leg2_steps = []
adjusted_body_move_central_steps = []
for i in range(len(body_move_leg1_steps)):
    leg1_upper, leg1_knee, leg1_ancle = body_move_leg1_steps[i]
    adjusted_leg1_knee = deg_curve(leg1_knee, leg1_ancle, last_deg2)
    adjusted_leg1_upper = deg_curve(leg1_upper, leg1_ancle, last_deg2)

    leg2_upper, leg2_knee, leg2_ancle = body_move_leg2_steps[i]
    adjusted_leg2_knee = deg_curve(leg2_knee, leg2_ancle, last_deg2)
    adjusted_leg2_upper = deg_curve(leg2_upper, leg2_ancle, last_deg2)

    adjusted_body_move_leg1_step = (adjusted_leg1_upper, adjusted_leg1_knee, leg1_ancle)
    adjusted_body_move_leg2_step = (adjusted_leg2_upper, adjusted_leg2_knee, leg2_ancle)

    central_position = ((adjusted_leg1_upper[0]+adjusted_leg2_upper[0])/2, adjusted_leg1_upper[1], adjusted_leg1_upper[2])

    adjusted_body_move_leg1_steps.append(adjusted_body_move_leg1_step)
    adjusted_body_move_leg2_steps.append(adjusted_body_move_leg2_step)
    adjusted_body_move_central_steps.append(central_position)

    last_deg2 = last_deg2-one_step_unit
    #print(last_deg2)

#print('adjusted_leg1:', adjusted_body_move_leg1_steps)
#print('adjusted_leg2:', adjusted_body_move_leg2_steps)
#print('adjusted_central:', adjusted_body_move_central_steps)
#print('original_leg1:', )

#print(deg_curve(adjusted_body_move_leg1_steps[-1][0], adjusted_body_move_leg1_steps[-1][2], 0))
#print(adjusted_body_move_leg1_steps[-1])

for_steps_leg2_upper = deg_curve(adjusted_body_move_leg2_steps[-1][0], adjusted_body_move_leg2_steps[-1][2], 0)
for_steps_leg2_ancle = deg_curve(adjusted_body_move_leg2_steps[-1][2], adjusted_body_move_leg2_steps[-1][2], 0)

#print('upper:', for_steps_leg2_upper)
#print('knee:', for_steps_leg2_ancle)

file.write('\n')
file.write('leg1 steps:{}\n'.format(leg1_step))
file.write('leg1 weight center movements:{}\n'.format(leg1_weight_center_movements))
file.write('still leg coordinates:{}\n'.format([still_leg_upper, still_leg_knee, still_leg_ancle]))
file.write('body particle:{}\n'.format(body_weight_center))
file.write('target position:{}\n'.format(target_position))
file.write('side move degrees:{}\n'.format(degrees))
file.write('adjusted leg1:{}\n'.format(adjusted_body_move_leg1_steps))
file.write('adjusted leg2:{}\n'.format(adjusted_body_move_leg2_steps))
file.write('adjusted central:{}\n'.format(adjusted_body_move_central_steps))

leg2_step, leg2_degs = leg_step_generator(how_far=body_forward+body_forward, time=one_action_time, body_fix=for_steps_leg2_upper, foot_origin=for_steps_leg2_ancle, lower_leg=below_leg, upper_leg=upper_leg)

status = 54
guage()
#print(leg1_degs)
first_upper_deg_in_steps, first_knee_deg_in_steps = leg2_degs[0][0], leg2_degs[0][1]
latest_deg_upper, latest_deg_knee = first_upper_deg_in_steps, first_knee_deg_in_steps
deg_changes_upper, deg_changes_knee = [0.0], [0.0]
for i in leg2_degs:
    change_upper = i[upper]-latest_deg_upper
    deg_changes_upper.append(change_upper)
    latest_deg_upper = i[upper]

    change_knee = i[knee]-latest_deg_knee
    deg_changes_knee.append(change_knee)
    latest_deg_knee = i[knee]

#print(deg_changes_upper)
#print(deg_changes_knee)

leg2_weight_center_movements = []
for i in range(len(leg2_step)):
    leg2_weight_center, anclezy, kneezy = leg_center(leg_num='A', ancle=leg2_step[i][2], knee=leg2_step[i][1], upper=leg2_step[i][0], ancle_mass=leg2_ancle_mass, knee_mass=leg2_knee_mass, upper_mass=leg2_upper_mass, degree_upper=deg_changes_upper[i], degree_knee=deg_changes_knee[i], upper_time=sleep_time, knee_time=sleep_time)
    manipulated = (leg2_step[i][0][0], leg2_weight_center[0][1], leg2_weight_center[0][0])
    leg2_weight_center_movements.append(manipulated)

still_leg_upper = (leg1_upper_particle[0], leg1_knee_particle[1], leg1_ancle_particle[2]+body_forward)
still_leg_knee = body_move_leg1_steps[-1][1]
still_leg_ancle = (leg1_ancle_particle[0], leg1_ancle_particle[1], leg1_ancle_particle[2]+body_forward)
still_leg_weight_center = weight_center_3D([still_leg_upper, still_leg_knee, still_leg_ancle], [leg1_upper_mass, leg1_knee_mass, leg1_ancle_mass])

body_weight_center = (body_weight_particle[0], body_weight_particle[1], body_weight_particle[2]+body_forward)
target_position = (leg1_ancle_particle[0], leg1_ancle_particle[1], leg1_ancle_particle[2]+body_forward)

degrees2 = [stablized_degree]   #body move final result!

right_leg_movement = []
left_leg_movement = []
central_movement = []
bowls = []
for a, i in enumerate(leg2_weight_center_movements):
    degree = 0
    #print('target')
    #print('for-degree', degree)
    match = 0
    bowl = []
    count = 5

    #print(count)
    while match == 0 and len(bowl) != count:
        last_deg = 0
        try:
            if degrees2[-1] != 0:
                last_deg = degrees2[-1]
        except:
            st = 'no last deg'
            #print('no last deg')
        adjust_left_leg = deg_curve(i, leg2_step[a][2], last_deg)
        adjust_right_leg = deg_curve(still_leg_weight_center, still_leg_ancle, last_deg)
        adjust_central = deg_curve(body_weight_center, ((leg2_step[a][2][0]+still_leg_ancle[0])/2, leg2_step[a][2][1]), last_deg)
        weight_center_changes, right_leg_change, left_leg_change, central_change = side_moves(right_leg=adjust_right_leg, left_leg=adjust_left_leg, left_foot=leg2_step[a][2], right_foot=still_leg_ancle, right_mass=leg1_upper_mass+leg1_knee_mass+leg1_ancle_mass, left_mass=leg2_upper_mass+leg2_knee_mass+leg2_ancle_mass, central_mass=body_weight_center_mass, central_particle=adjust_central, ancle_theta=degree*0.1, first_ancle_theta=last_deg, time=sleep_time*(a+1), wait_time=0)
        #print(weight_center_3D([i, still_leg_weight_center, body_weight_center], [5, 5, 3])[2])
        #print(last_deg)
        #print('degree:', degree*0.1+last_deg)
        '''if len(bowl)==0:
            target = 0
        else:
            target = bowl[-1]'''
        #print(target-degree*0.1, ',', target)
        #if abs(target-degree*0.1)>1:
        if abs(weight_center_changes[-1][0] - target_position[0]) < gap_allow*6:
            bowl.append(degree * 0.1 + last_deg)
            right_leg_movement.append(right_leg_change)
            left_leg_movement.append(left_leg_movement)
            central_movement.append(central_movement)

            if len(bowl) == count:
                match = 1

            if degree * 0.1 > 90:
                break

        degree -= 1


    #print(bowl)
    to_compare = degrees2[0]
    bowl2 = []
    for b in bowl:
        thing = abs(to_compare-b)
        bowl2.append(thing)
    minimum_one = min(bowl2)
    min_index = bowl2.index(minimum_one)
    the_degree = bowl[min_index]
    degrees2.append(the_degree)
    bowls.append(bowl)
    if (a+1)==(len(leg2_weight_center_movements)-1)/2:
        break

#print(degrees2)
#print('degree2:', len(degrees2))
#print(len(leg2_weight_center_movements))

for i in range(len(degrees2)-1):
    degrees2.append(degrees2[len(degrees2)-2-2*i])
#print('count', len(degrees2))
#print(bowls)
#print(degrees2)
latest_leg1_upper, latest_leg1_knee, latest_leg1_ancle = body_move_leg1_steps[-1]
latest_leg2_upper, latest_leg2_knee, latest_leg2_ancle = leg2_step[-1]
body_move_leg1_steps = [body_move_leg1_steps[-1]]

for i in range(50):
    ancle1 = latest_leg1_ancle
    upper1 = (latest_leg1_upper[0], latest_leg1_upper[1], latest_leg1_upper[2]+((i+1)*one_body_forward_unit))
    knee_points = get_knee_points(upper1, ancle1, upper_leg=upper_leg, lower_leg=below_leg)
    correct_point = 0
    for a in knee_points:
        degree, side = degree_finder(upper1, a, ancle1)
        if side == 'left':
            correct_point = a
    to_return = [upper1, correct_point, ancle1]
    body_move_leg1_steps.append(to_return)
    body_forward_current += one_body_forward_unit

body_forward_current2 = 0
body_move_leg2_steps = [(latest_leg2_upper, latest_leg2_knee, latest_leg2_ancle)]
for i in range(50):
    ancle1 = latest_leg2_ancle
    upper1 = (latest_leg2_upper[0], latest_leg2_upper[1], latest_leg2_upper[2]+((i+1)*one_body_forward_unit))
    knee_points = get_knee_points(upper1, ancle1, upper_leg=upper_leg, lower_leg=below_leg)
    correct_point = 0
    for a in knee_points:
        degree, side = degree_finder(upper1, a, ancle1)
        if side == 'left':
            correct_point = a
    to_return = [upper1, correct_point, ancle1]
    body_move_leg2_steps.append(to_return)
    body_forward_current2 += one_body_forward_unit

#print('body move 1:', body_move_leg1_steps)
#print('body move 2:', body_move_leg2_steps)
#print(len(body_move_leg1_steps), len(body_move_leg2_steps))

stablized_degree = 0
stablize_degree = 0
match_for_stablize = 0
instant_target = body_move_leg2_steps[-1][2]
leg1_weight_center_instant = weight_center_3D(body_move_leg1_steps[-1], [leg1_upper_mass, leg1_knee_mass, leg1_ancle_mass])
leg2_weight_center_instant = weight_center_3D(body_move_leg2_steps[-1], [leg2_upper_mass, leg2_knee_mass, leg2_ancle_mass])
central_position_edited = (body_weight_particle[0], body_weight_particle[1], body_weight_particle[2]+body_forward+body_forward)
while match_for_stablize == 0:
    instant_weight_center = weight_center_3D([leg1_weight_center_instant, leg2_weight_center_instant, central_position_edited], [leg1_upper_mass+leg1_knee_mass+leg1_ancle_mass, leg2_upper_mass+leg2_knee_mass+leg2_ancle_mass, body_weight_center_mass])
    adjusted_instant_weight_center = deg_curve(instant_weight_center, ((body_move_leg1_steps[-1][2][0]+body_move_leg2_steps[-1][2][0])/2, body_move_leg1_steps[-1][2][1]), stablize_degree*0.1)
    if abs(adjusted_instant_weight_center[0]-instant_target[0]) < gap_allow:
        stablized_degree = stablize_degree*0.1
        match_for_stablize = 1

    stablize_degree += 1

#print(stablized_degree)

one_step_unit = ((abs(stablized_degree)+abs(degrees2[-1])))/(len(body_move_leg1_steps))
last_deg2 = degrees2[-1]
adjusted_body_move_leg1_steps2 = []
adjusted_body_move_leg2_steps2 = []
adjusted_body_move_central_steps2 = []
for i in range(len(body_move_leg1_steps)):
    leg1_upper, leg1_knee, leg1_ancle = body_move_leg1_steps[i]
    adjusted_leg1_knee = deg_curve(leg1_knee, leg1_ancle, last_deg2)
    adjusted_leg1_upper = deg_curve(leg1_upper, leg1_ancle, last_deg2)

    leg2_upper, leg2_knee, leg2_ancle = body_move_leg2_steps[i]
    adjusted_leg2_knee = deg_curve(leg2_knee, leg2_ancle, last_deg2)
    adjusted_leg2_upper = deg_curve(leg2_upper, leg2_ancle, last_deg2)

    adjusted_body_move_leg1_step = (adjusted_leg1_upper, adjusted_leg1_knee, leg1_ancle)
    adjusted_body_move_leg2_step = (adjusted_leg2_upper, adjusted_leg2_knee, leg2_ancle)

    central_position = ((adjusted_leg1_upper[0]+adjusted_leg2_upper[0])/2, adjusted_leg1_upper[1], adjusted_leg1_upper[2])

    adjusted_body_move_leg1_steps2.append(adjusted_body_move_leg1_step)
    adjusted_body_move_leg2_steps2.append(adjusted_body_move_leg2_step)
    adjusted_body_move_central_steps2.append(central_position)

    last_deg2 = last_deg2+one_step_unit
    #print(last_deg2)

status = 78
guage()

#print('adjusted_leg1:', adjusted_body_move_leg1_steps2)
#print('adjusted_leg2:', adjusted_body_move_leg2_steps2)
#print('adjusted_central:', adjusted_body_move_central_steps2)
#print('original_leg1:', )

#print(deg_curve(adjusted_body_move_leg1_steps2[-1][0], adjusted_body_move_leg1_steps2[-1][2], 0))
#print(adjusted_body_move_leg1_steps[-1])

for_steps_leg1_upper = deg_curve(adjusted_body_move_leg1_steps2[-1][0], adjusted_body_move_leg1_steps2[-1][2], 0)
for_steps_leg1_ancle = deg_curve(adjusted_body_move_leg1_steps2[-1][2], adjusted_body_move_leg1_steps2[-1][2], 0)

#print('upper:', for_steps_leg1_upper)
#print('knee:', for_steps_leg1_ancle)

file.write('\n')
file.write('leg2 steps:{}\n'.format(leg2_step))
file.write('leg2 weight center movements:{}\n'.format(leg2_weight_center_movements))
file.write('still leg coordinates:{}\n'.format([still_leg_upper, still_leg_knee, still_leg_ancle]))
file.write('body particle:{}\n'.format(body_weight_center))
file.write('target position:{}\n'.format(target_position))
file.write('side move degrees:{}\n'.format(degrees2))
file.write('adjusted leg1:{}\n'.format(adjusted_body_move_leg1_steps))
file.write('adjusted leg2:{}\n'.format(adjusted_body_move_leg2_steps))
file.write('adjusted central:{}\n'.format(adjusted_body_move_central_steps))

leg3_step, leg3_degs = leg_step_generator(how_far=body_forward+body_forward, time=one_action_time, body_fix=for_steps_leg1_upper, foot_origin=for_steps_leg1_ancle, lower_leg=below_leg, upper_leg=upper_leg)

status = 90
guage()
#print(leg1_degs)
first_upper_deg_in_steps, first_knee_deg_in_steps = leg3_degs[0][0], leg3_degs[0][1]
latest_deg_upper, latest_deg_knee = first_upper_deg_in_steps, first_knee_deg_in_steps
deg_changes_upper, deg_changes_knee = [0.0], [0.0]
for i in leg3_degs:
    change_upper = i[upper]-latest_deg_upper
    deg_changes_upper.append(change_upper)
    latest_deg_upper = i[upper]

    change_knee = i[knee]-latest_deg_knee
    deg_changes_knee.append(change_knee)
    latest_deg_knee = i[knee]

#print(deg_changes_upper)
#print(deg_changes_knee)

leg1_weight_center_movements = []
for i in range(len(leg3_step)):
    leg2_weight_center, anclezy, kneezy = leg_center(leg_num='A', ancle=leg3_step[i][2], knee=leg3_step[i][1], upper=leg3_step[i][0], ancle_mass=leg1_ancle_mass, knee_mass=leg1_knee_mass, upper_mass=leg1_upper_mass, degree_upper=deg_changes_upper[i], degree_knee=deg_changes_knee[i], upper_time=sleep_time, knee_time=sleep_time)
    manipulated = (leg3_step[i][0][0], leg2_weight_center[0][1], leg2_weight_center[0][0])
    leg1_weight_center_movements.append(manipulated)

still_leg_upper = body_move_leg2_steps[-1][0]
still_leg_knee = body_move_leg2_steps[-1][1]
still_leg_ancle = body_move_leg2_steps[-1][2]
still_leg_weight_center = weight_center_3D([still_leg_upper, still_leg_knee, still_leg_ancle], [leg2_upper_mass, leg2_knee_mass, leg2_ancle_mass])

body_weight_center = (body_weight_particle[0], body_weight_particle[1], body_weight_particle[2]+body_forward+body_forward)
target_position = (leg2_ancle_particle[0], leg2_ancle_particle[1], leg2_ancle_particle[2]+body_forward+body_forward)

degrees3 = [stablized_degree]   #body move final result!

right_leg_movement = []
left_leg_movement = []
central_movement = []
bowls = []
for a, i in enumerate(leg1_weight_center_movements):
    degree = 0
    #print('for-degree', degree)
    match = 0
    bowl = []
    count = 10

    #print(count)
    while match == 0 and len(bowl) != count:
        last_deg = 0
        try:
            if degrees3[-1] != 0:
                last_deg = degrees3[-1]
        except:
            st = 'no last deg'
            #print('no last deg')
        adjust_right_leg = deg_curve(i, leg3_step[a][2], last_deg)
        adjust_left_leg = deg_curve(still_leg_weight_center, still_leg_ancle, last_deg)
        adjust_central = deg_curve(body_weight_center, ((leg3_step[a][2][0]+still_leg_ancle[0])/2, leg3_step[a][2][1]), last_deg)
        weight_center_changes, right_leg_change, left_leg_change, central_change = side_moves(right_leg=adjust_right_leg, left_leg=adjust_left_leg, right_foot=leg3_step[a][2], left_foot=still_leg_ancle, right_mass=leg1_upper_mass+leg1_knee_mass+leg1_ancle_mass, left_mass=leg2_upper_mass+leg2_knee_mass+leg2_ancle_mass, central_mass=body_weight_center_mass, central_particle=adjust_central, ancle_theta=degree*0.1, first_ancle_theta=last_deg, time=sleep_time*(a+1), wait_time=0)
        #print(weight_center_3D([i, still_leg_weight_center, body_weight_center], [5, 5, 3])[2])
        #print(last_deg)
        #print('degree:', degree*0.1+last_deg)
        '''if len(bowl)==0:
            target = 0
        else:
            target = bowl[-1]'''
        #print(target-degree*0.1, ',', target)
        #if abs(target-degree*0.1)>1:
        if abs(weight_center_changes[-1][0] - target_position[0]) < gap_allow*6:
            bowl.append(degree * 0.1 + last_deg)
            right_leg_movement.append(right_leg_change)
            left_leg_movement.append(left_leg_movement)
            central_movement.append(central_movement)

            if len(bowl) == count:
                match = 1

            if degree * 0.1 > 90:
                break

        degree += 1


    #print(bowl)
    to_compare = degrees2[0]
    bowl2 = []
    for b in bowl:
        thing = abs(to_compare-b)
        bowl2.append(thing)
    minimum_one = min(bowl2)
    min_index = bowl2.index(minimum_one)
    the_degree = bowl[min_index]
    degrees3.append(the_degree)
    bowls.append(bowl)
    if (a+1)==(len(leg1_weight_center_movements)-1)/2:
        break

status = 93
guage()
#print(degrees3)
#print('degree3:', len(degrees3))
#print(len(leg1_weight_center_movements))

for i in range(len(degrees3)-1):
    degrees3.append(degrees3[len(degrees3)-2-2*i])
#print('count', len(degrees3))
#print(bowls)
#print(degrees3)
#print(still_leg_upper, still_leg_knee, still_leg_ancle)

file.write('\n')
file.write('leg3 steps:{}\n'.format(leg3_step))
file.write('leg3 weight center movements:{}\n'.format(leg1_weight_center_movements))
file.write('still leg coordinates:{}\n'.format([still_leg_upper, still_leg_knee, still_leg_ancle]))
file.write('body particle:{}\n'.format(body_weight_center))
file.write('target position:{}\n'.format(target_position))
file.write('side move degrees3:{}\n'.format(degrees3))
file.write('adjusted leg1:{}\n'.format(adjusted_body_move_leg1_steps))
file.write('adjusted leg2:{}\n'.format(adjusted_body_move_leg2_steps))
file.write('adjusted central:{}\n'.format(adjusted_body_move_central_steps))

status = 100
guage()

print('\nTime Spent:{} Seconds'.format(time.time()-starttime))

file.close()
while True:
    answer = input('\nEXIT?:')
    if answer == 'yes':
        break
