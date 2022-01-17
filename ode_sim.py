import numpy as np
import matplotlib.pyplot as plt

def update_system(x, t, f, dt, sign = 1):
    x_dot = f(x)
    return np.add(x, np.multiply(x_dot, dt * sign)) , np.add(t, sign * dt) 

def simulate_system(x0, t0, f, dt, t_min, t_max):
    xs = [x0]
    ts = [t0]
    t = t0
    x = x0
    
    while t < t_max and np.abs(np.amax(x)) < 10**12:
        temp_x, temp_t = update_system(x, t, f, dt, sign = 1)
        # can be cleaned up by calling xs[-1], ts[-1] for updating
        xs.append(temp_x)
        ts.append(temp_t)
        t = temp_t
        x = temp_x
        
    t = t0
    x = x0
    
    while t > t_min and np.abs(np.amax(x)) < 10**12:
        temp_x, temp_t = update_system(x, t, f, dt, sign = -1)
        # can be cleaned up by calling xs[0], ts[0] for updating
        xs.insert(0, temp_x)
        ts.insert(0, temp_t)
        t = temp_t
        x = temp_x
    
    return ts, xs

def draw_phase_curves(x0s, t0s, f, 
                      dt, t_min, t_max, name
                      xrange = [-3, 3], yrange = [-1.5, 1.5],
                      type = 0, title = "", size = 5):
    # generally used to visualize single graphs - subplots generated separately
    plt.figure(figsize = (size * (yrange[1] - yrange[0]) / (xrange[1] - xrange[0]), size))
    for i in range(len(x0s)):
        times, positions = simulate_system(x0s[i], t0s[i], f, dt, t_min, t_max)
        if type == 0:
            plt.plot(times, positions, color = "orange")
        elif type == 1:
            #plt.plot(positions[0], positions[1], color = "orange")
            #print(positions)
            trans_pos = np.transpose(positions)
            plt.plot(trans_pos[0], trans_pos[1], color = "orange")
        else:
            print("draw_phase_curves: unsupported type")
            return

    plt.title(title)
    plt.ylim(yrange[0], yrange[1])
    plt.xlim(xrange[0], xrange[1])
    plt.axvline(x = 0, c = "black")
    plt.axhline(y = 0, c = "black")
    plt.savefig(name)
    
def lotka_volterra(k, l, a, b):
    return lambda x: [k * x[0] - a * x[0] * x[1], - l * x[1] + b * x[0] * x[1]]

n = 5
lv_init_pos = np.transpose([np.linspace(2, 5, n), np.ones(n)])
draw_phase_curves(np.transpose([np.linspace(2, 10, n), np.ones(n)]), np.zeros(n), 
                  lotka_volterra(6, 4, 2, 2), 10**(-5), 0, 10, name = "l_v.png",
                  xrange = [0, 15], yrange = [0, 15], type = 1, title = "lv sim", size = 10)
