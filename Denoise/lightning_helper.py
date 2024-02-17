from abc import abstractmethod
import pytorch_lightning as pl

class LightningNetwork(pl.LightningModule):
    def __init__(self, network, optimizer):
        super(LightningNetwork, self).__init__()
        self.network = network
        self.optimizer = optimizer

    @abstractmethod
    def training_step(self, batch, batch_idx):
        pass

    def configure_optimizers(self):
        return self.optimizer
