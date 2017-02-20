/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 * @ Anna
 * This class contains the pre-computes gradients of the volume. This means calculates the gradient
 * at all voxel positions, and provides functions
 * to get the gradient at any position in the volume also continuous..
*/
public class GradientVolume {

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
        
    }

    public VoxelGradient getGradient(int x, int y, int z) {
         return data[x + dimX * (y + dimY * z)];
    }

    private VoxelGradient interpolate(VoxelGradient g0, VoxelGradient g1, float factor) {
        /* Implemented: this function linearly interpolates gradient vector g0 and g1 given the factor (t) 
            the resut is given at result. You can use it to tri-linearly interpolate the gradient */
        VoxelGradient result = new VoxelGradient(g0.x*(1-factor) + g1.x*factor, g0.y*(1-factor) + g1.y*factor,g0.z*(1-factor) + g1.z*factor);
        return result; 
    }
    
    public VoxelGradient getGradientNN(double[] coord) {
        /* Nearest neighbour interpolation applied to provide the gradient */
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }

        int x = (int) Math.round(coord[0]);
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
        
        return getGradient(x, y, z);
    }

    
    public VoxelGradient getGradient(double[] coord) {
    /* Implemented: Returns trilinear interpolated gradient based on the precomputed gradients. 
     *   Use function interpolate. Use getGradientNN as bases 
    */
   
    //Get gradient of the interpolated values
     if (coord[0] < 0 || coord[0] > volume.getDimX() - 1 || coord[1] < 0 || coord[1] > volume.getDimY() - 1
                || coord[2] < 0 || coord[2] > volume.getDimZ() - 1) {
            return getGradient(0,0,0);
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        if (x > volume.getDimX() - 2
                || y > volume.getDimY() - 2
                || z > volume.getDimZ() - 2) {
            return getGradient(x,y,z);
        }
        VoxelGradient val0 = getGradient(x, y, z);
        VoxelGradient val1 = getGradient(x + 1, y, z);
        
        VoxelGradient val2 = getGradient(x, y + 1, z);
        VoxelGradient val3 = getGradient(x + 1, y + 1, z);
        
        VoxelGradient val4 = getGradient(x, y, z + 1);
        VoxelGradient val5 = getGradient(x + 1, y, z + 1);
        
        VoxelGradient val6 = getGradient(x, y + 1, z + 1);
        VoxelGradient val7 = getGradient(x + 1, y + 1, z + 1);
        
        float alpha = (float) (coord[0] - x);
        float beta = (float) (coord[1] - y);
        float gamma = (float) (coord[2] - z);
        
        // four linear interpolation
        VoxelGradient val01 = interpolate(val0,val1,alpha);
        VoxelGradient val23 = interpolate(val2,val3,alpha);
        VoxelGradient val45 = interpolate(val4,val5,alpha);
        VoxelGradient val67 = interpolate(val6,val7,alpha);
              
        //two bilinear interpolations
        VoxelGradient val0123 = interpolate(val01,val23,beta);
        VoxelGradient val4567 = interpolate(val45,val67,beta);
         
        
        //tri linear interpolation
        VoxelGradient finalVal = interpolate(val0123,val4567,gamma);
         

        return finalVal;
  
    }
    
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }

    public VoxelGradient getVoxel(int i) {
        return data[i];
    }

    public int getDimX() {
        return dimX;
    }

    public int getDimY() {
        return dimY;
    }

    public int getDimZ() {
        return dimZ;
    }

    private void compute() {
        /*implemented: compute the gradient of contained in the volume attribute */
        for (int x = 0; x < volume.getDimX(); x ++) {
            for (int y = 0; y < volume.getDimY(); y++) {
                for (int z = 0; z < volume.getDimZ(); z++) {
                    VoxelGradient gradient;
                    if (x == 0 || x == volume.getDimX() - 1 
                            || y == 0 || y == volume.getDimY() - 1
                            || z == 0 || z == volume.getDimZ() - 1) {
                        gradient = new VoxelGradient(0, 0, 0);
                    } else {
                        gradient = new VoxelGradient(
                               0.5f * (volume.getVoxel(x + 1, y, z) - volume.getVoxel(x - 1, y, z))
                               , 0.5f * (volume.getVoxel(x, y + 1, z) - volume.getVoxel(x, y - 1, z))
                               , 0.5f * (volume.getVoxel(x, y, z + 1) - volume.getVoxel(x, y, z - 1))
                       );
                    }
                    this.setGradient(x, y, z, gradient);
                }
            }
        }  
     
    }
    
    public double getMaxGradientMagnitude() {
        /* implemented: Returns the maximum gradient magnitude*/
        if(maxmag >= 0)
        { return maxmag;}
        else {
        double magnitude = data[0].mag;
        for(int i=0; i<data.length; i++)
        { magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;}
        maxmag = magnitude;
        return magnitude;
        }
    }
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
}