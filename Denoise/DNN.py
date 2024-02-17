import torch
import torch.nn as nn
import torch.nn.functional as F

class TorchMLP(nn.Module):
    def __init__(self, dimensions, bias=True, activation_function=nn.Tanh(), bn=None, dropout=None, final_nonlinearity=False):
        super(TorchMLP, self).__init__()
        # setup ModuleList
        self.model  = nn.ModuleList()
        for i in range(len(dimensions)-1):
            self.model.append(nn.Linear(dimensions[i], dimensions[i+1], bias=bias))
            if bn is not None and i!=len(dimensions)-2:
                self.model.append(bn(dimensions[i+1]))
                if dropout is not None:
                    self.model.append(dropout)
            if i==len(dimensions)-2: break
            self.model.append(activation_function)
        if final_nonlinearity:
            if bn is not None:
                self.model.append(bn(dimensions[-1]))
            self.model.append(activation_function)
        # weight init
        self.model.apply(self.xavier_init)

    def xavier_init(self, m):
        if type(m) == nn.Linear:
            torch.nn.init.xavier_uniform_(m.weight)
            m.bias.data.fill_(0.0)

    def forward(self, x):
        for i, l in enumerate(self.model): 
            x = l(x)
        return x

class Sine(nn.Module):
    def __init__(self, ):
        super(Sine, self).__init__()
    def forward(self, x):
        return torch.sin(x)
