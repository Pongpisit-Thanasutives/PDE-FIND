from DNN import *

def count_parameters(torch_model, onlyif_requires_grad=True):
    if onlyif_requires_grad:
        return sum(p.numel() for p in torch_model.parameters() if p.requires_grad)
    return sum(p.numel() for p in torch_model.parameters())

class PhysicalConstraintCalculator(nn.Module):
    def __init__(self, symbolic_module, basic_vars, init_coefficients=None, learnable_coefficients=False):
        super(PhysicalConstraintCalculator, self).__init__()
        self.symbolic_module = symbolic_module
        self.basic_vars = basic_vars

        self.coefficients = init_coefficients
        self.learnable_coefficients = learnable_coefficients

        if self.coefficients is None:
            self.coefficients = torch.ones(len(symbolic_module.sympy())).float()
        else:
            self.coefficients = torch.tensor(data=self.coefficients).float()
        self.coefficients = nn.Parameter(self.coefficients).requires_grad_(self.learnable_coefficients)

        # printing
        if self.learnable_coefficients: print("Learnable coefficients:", self.coefficients)
        else: print("NOT learnable coefficients:", self.coefficients)
        print(symbolic_module.sympy())
        print("Basic variables:", self.basic_vars)

    def set_learnable_coefficients(self, learn):
        self.coefficients.requires_grad_(learn)

    def forward(self, input_dict):
        return self.symbolic_module(**input_dict)

class PINN(nn.Module):
    def __init__(self, solver, physics_calculator, lb, ub, 
                 domain_dimension=None, weak_pde_lib=None, effective_indices=None):
        super(PINN, self).__init__()
        self.solver = solver
        self.physics_calculator = physics_calculator
        self.lb = lb
        self.ub = ub
        # Only to use weak_loss
        # spatial x temporal
        self.domain_dimension = domain_dimension
        self.weak_pde_lib = weak_pde_lib
        self.effective_indices = effective_indices
        self.weak_coeff_buffer = None
        
    def forward(self, x, t):
        return self.solver(self.input_normalize(torch.cat([x, t],  dim=-1)))

    def calculate_physics(self, x, t):
        u = self.forward(x, t)
        u_t = self.gradients(u, t)[0]
        u_1 = self.gradients(u, x)[0]
        u_11 = self.gradients(u_1, x)[0]
        physics = self.physics_calculator({nameof(u):u, 
                                           nameof(u_1):u_1, 
                                           nameof(u_11):u_11})
        
        return u, u_t, physics
    
    def loss(self, x, t, y_input):
        u, u_t, physics = self.calculate_physics(x, t)
        coeff = self.physics_calculator.coefficients
        physics = (physics*coeff).sum(axis=-1)
        mse = F.mse_loss(u, y_input, reduction='mean')
        l_eq = F.mse_loss(u_t, physics, reduction='mean')
        return torch.add(mse, l_eq)
    
    def weak_loss(self, x, t, y_input):
        u, u_t, physics = self.calculate_physics(x, t)
        coeff = torch.tensor(self.weak_coefficients(u)).float()
        physics = (physics*coeff).sum(axis=-1)
        mse = F.mse_loss(u, y_input, reduction='mean')
        l_eq = F.mse_loss(u_t, physics, reduction='mean')
        return torch.add(mse, l_eq)
    
    def weak_form(self, u):
        pred = u.reshape(self.domain_dimension[1], 
                         self.domain_dimension[0]).T.detach().numpy()
        pred = np.expand_dims(pred,-1)
        X_weak = self.weak_pde_lib.fit_transform(pred)
        y_weak = self.weak_pde_lib.convert_u_dot_integral(pred)
        return X_weak, y_weak
    
    def weak_coefficients(self, u):
        np.random.seed(0)
        X_weak, y_weak = self.weak_form(u)
        X_weak = X_weak[:, self.effective_indices]
        self.weak_coeff_buffer = np.linalg.lstsq(X_weak, y_weak, rcond=None)[0].flatten()
        return self.weak_coeff_buffer
    
    def gradients(self, func, x):
        return grad(func, x, create_graph=True, retain_graph=True, 
                    grad_outputs=torch.ones(func.shape))

    def input_normalize(self, inp):
        return -1.0+2.0*(inp-self.lb)/(self.ub-self.lb)
