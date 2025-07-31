import numpy as np
import csv

# 1D Finite Volume Method for shallow water equations
def main():
    # Physical parameters
    g = 9.81  # gravity (m/s^2)
    L = 14000.0  # domain length (m)
    final_time = 9117.5  # simulation end time (s)

    # Discretization parameters
    N = 500  # number of cells
    dx = L / N
    cfl = 0.9

    # Computational grid (cell centers)
    x = np.linspace(dx/2, L - dx/2, N)

    # Bathymetry and initial conditions functions
    def H(x):
        return 50.5 - 40 * x / L - 10 * np.sin(np.pi * (4 * x / L - 0.5))

    def zb(x):
        return H(0) - H(x)

    # Analytical solutions
    def h_analytical(x, t):
        return H(x) + 4 - 4 * np.sin(np.pi * (4 * t / 86400.0) + 0.5)

    def u_analytical(x, t):
        h_a = h_analytical(x, t)
        return (x - L) * np.pi / (5400 * h_a) * np.cos(np.pi * (4 * t / 86400.0 + 0.5))

    # Initialize variables
    h = H(x)
    u = np.zeros_like(x)
    zb_arr = zb(x)
    q = h * u  # discharge

    # Time step based on CFL
    max_speed = np.max(np.abs(u) + np.sqrt(g * h))
    dt = cfl * dx / (max_speed if max_speed > 0 else np.sqrt(g * np.max(h)))
    num_steps = int(np.ceil(final_time / dt))
    dt = final_time / num_steps  # adjust dt to hit final_time exactly


    # test = q**2 / h + 0.5 * g * h**2
    # print("Initial test value:", test)  # Debugging output
    # test = 0.5 * g * h**2
    # print(q**2)
    # Finite volume update (Lax-Friedrichs for illustration)
    for n in range(num_steps):
        t = n * dt
        if n % 100 == 0:  # Print every 100 steps
            print(f"Step {n}, Time {t:.2f}, dt {dt:.4f}, h: {h[:5]}, q: {q[:5]}")
        # fluxes
        Fh = q
        Fq = q**2 / h + 0.5 * g * h**2

        # Lax-Friedrichs update
        h_new = 0.5 * (np.roll(h, 1) + np.roll(h, -1)) - dt / (2 * dx) * (np.roll(Fh, -1) - np.roll(Fh, 1))
        q_new = 0.5 * (np.roll(q, 1) + np.roll(q, -1)) - dt / (2 * dx) * (np.roll(Fq, -1) - np.roll(Fq, 1))

        # Apply boundary conditions
        # h(0,t) = H(0) + 4 - 4 sin[pi(4t/86400)+1/2]
        h_bc_left = H(0) + 4 - 4 * np.sin(np.pi * (4 * (t+dt) / 86400.0) + 0.5)
        h_new[0] = h_bc_left
        q_new[ -1] = 0.0  # u(L,t) = 0 => q(L) = h(L)*0

        # Update
        h, q = h_new, q_new

    # Final velocity
    u = q / h

    # Analytical at final time
    h_a = h_analytical(x, final_time)
    u_a = u_analytical(x, final_time)
    zb_final = zb_arr

    # Write results to CSV
    with open('../results/results_fvm.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Header
        writer.writerow(['x_node', 'x_m', 'h_plus_zb', 'zb', 'h', 'u', 'h_analytical', 'u_analytical', 'h_analytical_plus_zb'])
        # Data rows
        for i in range(N):
            writer.writerow([
                i,
                x[i],
                h[i] + zb_final[i],
                zb_final[i],
                h[i],
                u[i],
                h_a[i],
                u_a[i],
                h_a[i] + zb_final[i]
            ])

    print("Simulation complete. Results written to ../results/results_fvm.csv")

if __name__ == '__main__':
    main()
