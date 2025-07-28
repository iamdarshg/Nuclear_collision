import numpy as np

# Re-run after kernel reset: regenerate and save PDF data
x_grid = np.logspace(-5, 0, 100)       # x from 1e-5 to 1
Q2_grid = np.logspace(-10, 3, 100)       # Q^2 from 1 to 1000 GeV^2

# Create mock PDF values for each parton flavor
def mock_pdf_func(flavor_index):
    return np.outer(np.exp(-flavor_index * np.log(x_grid + 1e-10)), 
                    1.0 / (1 + np.log(Q2_grid + 1)))

flavors = ['u', 'd', 's', 'c', 'b', 'g']
pdfs = {flavor: mock_pdf_func(i + 1) for i, flavor in enumerate(flavors)}

# Save everything to a .npz file
np.savez('./data/pdf_data.npz', x_grid=x_grid, Q2_grid=Q2_grid, **pdfs)

# This code generates a mock PDF dataset and saves it in a format that can be used by the APFEL library.