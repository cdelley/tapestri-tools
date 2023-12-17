import h5py as h5
import anndata as ad
import pandas as pd


class read_h5(object):
    """
    simple class to scan the content of a .h5 file
    based on: 
    https://stackoverflow.com/questions/43371438/how-to-inspect-h5-file-in-python/43374773#43374773
    
    """
    
    def __init__(self, path):
        self.nodes = []
        self.filename = path
    
    def _scan_node(self, g, tabs=0, recursive=True, tab_step=2, print_nodes=False) -> None:
        if print_nodes:
            print(' ' * tabs, g.name)
        for k, v in g.items():
            if isinstance(v, h5.Dataset):
                self.nodes.append(v.name)
                if print_nodes:
                    print(' ' * tabs + ' ' * tab_step + ' -', v.name)
            elif isinstance(v, h5.Group) and recursive:
                self._scan_node(v, tabs=tabs + tab_step, print_nodes=print_nodes)
    
    def scan_hdf5(self, recursive=True, tabs=0, tab_step=2, print_nodes=True) -> None:
        with h5.File(self.filename, 'r') as f:
            self._scan_node(f, tabs=tabs, recursive=recursive, tab_step=tab_step, print_nodes=print_nodes)

class read_tapestry(read_h5):
    """
    Tapestri outputs are .h5 file containers.
    """
    
    def __init__(self, path):
        super().__init__(path)
            
    def make_anndata(self, filtered=False) -> None:
        
        if filtered:
            datagroup = '/assays'
        else:
            datagroup = '/all_barcodes'
        
        with h5.File(self.filename, 'r') as f:
            #dset = f['/assays/dna_variants/ra/barcode'][()]
            #barcodes = dset.astype(str)

            var = {}
            for k in f['assays/dna_read_counts/ca'].keys():
                data = f['assays/dna_read_counts/ca/'+k][()]
                try:
                    data.astype(float)
                    var[k] = data.astype(data.dtype)
                except ValueError:
                    var[k] = data.astype(str)
            var['CHROM'] = var['CHROM'].astype(str)
            
            
            obs = {}
            for k in f[datagroup+'/dna_read_counts/ra'].keys():
                data = f[datagroup+'/dna_read_counts/ra/'+k][()]
                try:
                    data.astype(float)
                    obs[k] = data.astype(data.dtype)
                except ValueError:
                    obs[k] = data.astype(str)
            
            uns = {}
            for k in f[datagroup+'/dna_read_counts/metadata'].keys():
                data = f[datagroup+'/dna_read_counts/metadata/'+k][()]
                try:
                    Number = data.astype(float)
                    if Number.is_integer():
                        Number = int(Number)
                    uns[k] = Number
                except ValueError:
                    uns[k] = data.astype(str)
                    
            sel = {}
            for k in f['/assays/dna_variants/ra/'].keys():
                data = f['/assays/dna_variants/ra/'+k][()]
                try:
                    data.astype(float)
                    sel[k] = data.astype(data.dtype)
                except ValueError:
                    sel[k] = data.astype(str)
                    
            for k in f[datagroup+'/dna_read_counts/layers'].keys():
                data = f[datagroup+'/dna_read_counts/layers/'+k][()].astype(int)
        
        var=pd.DataFrame(var)
        obs=pd.DataFrame(obs)
        adata = ad.AnnData(
            X=data,
            var=var,
            obs=obs,
            uns=uns,
        )
        
        adata.obs_names = adata.obs['barcode']
        del adata.obs['barcode']
        adata.var_names = adata.var['id']
        del adata.var['id']
        
        self.adata = adata
        self.cell_bc = sel['barcode'] 

