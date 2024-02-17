from DNN import *

class DeepAutoencoder(nn.Module):
    def __init__(self, dimensions, bias=True, activation_function=nn.Sigmoid(), bn=None, dropout=None, final_nonlinearity=(False, False)):
        super(DeepAutoencoder, self).__init__()
        self.encoder = TorchMLP(dimensions=dimensions, bias=bias, activation_function=activation_function, bn=bn, dropout=dropout, final_nonlinearity=final_nonlinearity[0])
        self.decoder = TorchMLP(dimensions=dimensions[::-1], bias=bias, activation_function=activation_function, bn=bn, dropout=dropout, final_nonlinearity=final_nonlinearity[1])
    def forward(self, X):
        return self.decoder(self.encoder(X))
    def loss(self, X, loss_constant):
        recon = self.forward(X)
        return loss_constant*(torch.mean(torch.square(X-recon)))

class RobustL21Autoencoder(nn.Module):
    def __init__(self, dimensions, bias=True, activation_function=nn.Sigmoid(), bn=None, dropout=None, lambda_=1):
        self.ae = DeepAutoencoder(dimensions=dimensions, bias=bias, activation_function=activation_function, bn=bn, dropout=dropout)
        self.lambda_ = lambda_
