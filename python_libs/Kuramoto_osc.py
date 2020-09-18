from tqdm import tqdm
from time import time
import tensorflow as tf
import numpy as np
# import tensorflow.compat.v1 as tf
# tf.disable_v2_behavior()
print(tf.__version__)

class KuramotoOsc:
#     def __init__(self, N, adj_list, k, w0=None, theta0=None, dt=0.001):
    """
    #---------------descripts---------------#
    run Kuramoto oscillator with tensorflow 1.15.1
    #---------------arguments---------------#
    N: the # of oscillators, scalar
    k: connection strength
    adj_list: adjacency list of the network, dictionary
    w0 (optional): initial w
    theat0 (optional): initial theat
    """
    def __init__(self, **kwargs):
        self._set_constant(kwargs)
    
    def run(self, tmax, dt=0.001, desc="run"):
        # init tensor_obj
        self._init_run_obj(tmax, dt)
        # create tensor_obj
        self._init_tensor_obj()
        # run model
        with tf.Session() as sess:
            bar = tqdm(range(self.itr_all), desc=desc, ncols=100, mininterval=0.2)
            n_r = 0
            for i in bar:
                self.theta_list[i+1, np.newaxis, :] = sess.run(self.theta_next, feed_dict={self.theta_tf: self.theta_list[i, np.newaxis, :]})
                self.theta_list[i+1, :] = np.round(self.theta_list[i+1, :], 10)
                if (i+1) % 1000 == 0:
                    self.get_coh(self.theta_list[n_r:(i+2)])
                    n_r = i+2
            if i > n_r:
                self.get_coh(self.theta_list[n_r:])

    def get_coh(self, x):
        # calculate r and \psi
        temp = np.average(np.exp(1j * x), 1)
        r_temp = abs(temp)
        self.r = np.concatenate([self.r, abs(temp)])
        self.psi = np.concatenate([self.psi, np.log(temp/r_temp).imag])

    def _set_constant(self, kwargs):
        # check keys
        keys = list(kwargs.keys())
        key_need_sets = ["N", "k", "adj", "theta0", "w0"]
        for k in key_need_sets:
            if k in keys:
                keys.remove(k)
        if keys: # not empty
            print(keys, "another keys exits")
            return 0
        # init variables
        self.N = kwargs["N"]
        self.k = kwargs["k"]
        
        if "theta0" not in kwargs.keys():
            self.theta0 = np.random.normal(loc=0, scale=5, size=[1, self.N])
        else:
            self.theta0 = kwargs["theta0"]
        if "w0" not in kwargs.keys():
            self.w0 = abs(np.random.normal(loc=1, scale=0, size=[1, N]))
        else:
            self.w0 = kwargs["w0"]
        # check adjacency list or mat
        adj = kwargs["adj"]
        if type(adj) == dict:
            self.adj_list = adj
            # creaet adjacency matrix
            self.K_mat = np.zeros([self.N, self.N])
            for i in range(self.N):
                for j in self.adj_list[i]:
                    self.K_mat[i, j] = 1
            self.K_mat *= self.k
        else:
            self.adj_list = None
            self.K_mat = adj
            
    def _init_tensor_obj(self):
        # create tensorflow model
        self.theta_tf = tf.placeholder(tf.float32, shape=(1, self.N))
        self.K_tf = tf.constant(self.K_mat, dtype=tf.float32)
        self.w0_tf = tf.constant(self.w0, dtype=tf.float32)
        self.theta_next = self._update_obj()
        
    def _init_run_obj(self, tmax=10, dt=0.01):
        # time
        self.tmax = tmax
        self.dt = dt
        self.itr_all = int(self.tmax / self.dt)
        self.t = np.arange(0, self.tmax+self.dt/2, self.dt)
        # To save theat
        self.theta_list = np.zeros([self.itr_all+1, self.N]) # save theta
        self.theta_list[0, :] = self.theta0
        self.r = []
        self.psi = []
    
    def _update_obj(self):
        # use RK4 method
        dx1 = self.Kuramoto_eq(self.theta_tf)
        dx2 = self.Kuramoto_eq(self.theta_tf + 1/2*self.dt * dx1)
        dx3 = self.Kuramoto_eq(self.theta_tf + 1/2*self.dt * dx2)
        dx4 = self.Kuramoto_eq(self.theta_tf + self.dt * dx3)
        return self.theta_tf + 1/6*self.dt*(dx1 + 2*dx2 + 2*dx3 + dx4)
    
    def Kuramoto_eq(self, x):
        # x, [1, N] size
        dsin = tf.math.sin(tf.subtract(x, tf.transpose(x)))
        return tf.add(tf.reduce_sum(tf.multiply(self.K_tf, dsin), 1), self.w0_tf)

    def test_time_step(self):
        # postulate 1e-6 return the correct answer
        self.dt_sets = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 5e-5] 
        # self.dt_sets = [(0.1)**(i+1) for i in range(6)] # 0.1, 0.01, 0.001, 0.0001, 0.00001
        tmax = 1
        self.prec_x = []
        self.prec_r = []
        for dt in self.dt_sets:
            self.run(tmax, dt, desc="dt=%f"%(dt))
            self.prec_x.append(self.theta_list[-1, :]) # end points
            self.prec_r.append(self.r[-1])
        self.prec = []
        for i in range(len(self.dt_sets)-1):
            prec_dx = np.average(self.prec_x[i] - self.prec_x[-1])
            prec_dr = self.prec_r[i] - self.prec_r[-1]
            self.prec.append(prec_dx**2 + prec_dr**2)

if __name__ == "__main__":
    N = 2000
    adj_list = dict()
    edges_all = [i for i in range(N)]
    for i in range(N):
        edges = edges_all.copy()
        edges.pop(i)
        adj_list[i] = edges
    
    theta0 = np.random.normal(loc=0, scale=5, size=[1, N])
    w0 = abs(np.random.normal(loc=1, scale=0, size=[1, N]))

    osc = KuramotoOsc(N=N, adj_list=adj_list, k=1., theta0=theta0, w0=w0)
    osc.run(tmax=10, dt=0.001)

    plt.figure(dpi=200)
    for i in range(5):
        plt.plot(osc.t, np.sin(osc.theta_list[:, i]), "k", lw=0.5)
    plt.show()
