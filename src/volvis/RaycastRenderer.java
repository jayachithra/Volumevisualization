/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 * @Anna
 * This class has the main code that generates the raycasting result image. 
 * The connection with the interface is already given.  
 * The different modes mipMode, slicerMode, etc. are already correctly updated
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    private boolean mipMode = false;
    private boolean slicerMode = true;
    private boolean compositingMode = false;
    private boolean tf2dMode = false;
    private boolean shadingMode = false;
    
    /** Shading coefficients **/
    // ambient reflection coefficient, assuming light source is white
    TFColor SHADING_AMBIENT_COEFF = new TFColor(0.1, 0.1, 0.1, 1.0);
    // diffuse reflection coefficient
    double SHADING_DIFF_COEFF = 0.7;
    // specular reflection coefficient
    double SHADING_SPEC_COEFF = 0.2;
    // exponent used to approximate highligh
    double SHADING_ALPHA = 10;
    
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        tFunc.setTestFunc();
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     
    public void setShadingMode(boolean mode) {
        shadingMode = mode;
        changed();
    }
    
    public void setMIPMode() {
        setMode(false, true, false, false);
    }
    
    public void setSlicerMode() {
        setMode(true, false, false, false);
    }
    
    public void setCompositingMode() {
        setMode(false, false, true, false);
    }
    
    public void setTF2DMode() {
        setMode(false, false, false, true);
    }
    
    private void setMode(boolean slicer, boolean mip, boolean composite, boolean tf2d) {
        slicerMode = slicer;
        mipMode = mip;
        compositingMode = composite;
        tf2dMode = tf2d;        
        changed();
    }
    
        
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }
    
    private boolean intersectLinePlane(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection) {

        double[] tmp = new double[3];

        for (int i = 0; i < 3; i++) {
            tmp[i] = plane_pos[i] - line_pos[i];
        }

        double denom = VectorMath.dotproduct(line_dir, plane_normal);
        if (Math.abs(denom) < 1.0e-8) {
            return false;
        }

        double t = VectorMath.dotproduct(tmp, plane_normal) / denom;

        for (int i = 0; i < 3; i++) {
            intersection[i] = line_pos[i] + t * line_dir[i];
        }

        return true;
    }

    private boolean validIntersection(double[] intersection, double xb, double xe, double yb,
            double ye, double zb, double ze) {

        return (((xb - 0.5) <= intersection[0]) && (intersection[0] <= (xe + 0.5))
                && ((yb - 0.5) <= intersection[1]) && (intersection[1] <= (ye + 0.5))
                && ((zb - 0.5) <= intersection[2]) && (intersection[2] <= (ze + 0.5)));

    }

    private void intersectFace(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection,
            double[] entryPoint, double[] exitPoint) {

        boolean intersect = intersectLinePlane(plane_pos, plane_normal, line_pos, line_dir,
                intersection);
        if (intersect) {

            double xpos0 = 0;
            double xpos1 = volume.getDimX();
            double ypos0 = 0;
            double ypos1 = volume.getDimY();
            double zpos0 = 0;
            double zpos1 = volume.getDimZ();

            if (validIntersection(intersection, xpos0, xpos1, ypos0, ypos1,
                    zpos0, zpos1)) {
                if (VectorMath.dotproduct(line_dir, plane_normal) > 0) {
                    entryPoint[0] = intersection[0];
                    entryPoint[1] = intersection[1];
                    entryPoint[2] = intersection[2];
                } else {
                    exitPoint[0] = intersection[0];
                    exitPoint[1] = intersection[1];
                    exitPoint[2] = intersection[2];
                }
            }
        }
    }
    
    /**
     * Back-to-front compositing  
     */
    private int traceRayCompositingB2F(double[] entryPoint, double[]exitPoint, double[] viewVec,double sampleStep){
        double total_dis = VectorMath.distance(entryPoint, exitPoint);
        double entry_exit_vector[] = new double[3];
        double current_point[] = new double[3];
        short current_intensity = 0;
        TFColor basic_color = new TFColor(); //the basic color of each point from the TF
        TFColor aug_color = new TFColor(); // the augmented color after using compositing method
        aug_color.a = 1;
        aug_color.r = 0;
        aug_color.g = 0;
        aug_color.b = 0;
     //vector difference between exit Point and entry point define a plane through the origin for x,y and z plane
     //perpendicular to the entry_exit_vector
        VectorMath.setVector(entry_exit_vector, exitPoint [0]-entryPoint[0], exitPoint[1]-entryPoint[1], exitPoint[2]-entryPoint[2]);
        for (double current_dis = total_dis; current_dis >0; current_dis -= sampleStep){
            for(int i = 0; i<3; i++){
                current_point[i] = (current_dis/total_dis)*entry_exit_vector[i]+entryPoint[i];
            }
            current_intensity = volume.getVoxelInterpolate(current_point);
            basic_color = tFunc.getColor(current_intensity);
            float alpha = (float) (1 - Math.pow(1-basic_color.a, sampleStep)); // true opacity after the interpolation
            aug_color.r = basic_color.r * alpha +(1-alpha)*aug_color.r; // in Back-to-front compositing, the augmented color are composed with two parts:
            aug_color.g = basic_color.g *alpha + (1-alpha)*aug_color.g;  // 1. the basic color multiplied by its opacity just at the point(alpha); 
            aug_color.b = basic_color.b * alpha + (1-alpha)*aug_color.b;  // 2. former color transmitted through the point (1-alpha)
        }
        current_intensity = volume.getVoxelInterpolate(entryPoint);
        basic_color = tFunc.getColor(current_intensity);
        float alpha = (float) (1 - Math.pow(1-basic_color.a, sampleStep));
        aug_color.r = basic_color.r * alpha +(1-alpha)*aug_color.r;
        aug_color.g = basic_color.g *alpha + (1-alpha)*aug_color.g;
        aug_color.b = basic_color.b * alpha + (1-alpha)*aug_color.b;
       // aug_color.a = basic_color.a;

        return doublesToColor(1,aug_color.r,aug_color.g,aug_color.b);

    }

    /**
     * Front-to-back compositing 
    **/
   private int traceRayCompositingF2B(double[] entryPoint, double[]exitPoint, double[] viewVec,double sampleStep){
        double total_dis = VectorMath.distance(entryPoint, exitPoint);
       double entry_exit_vector[] = new double[3];
       double current_point[] = new double[3];
       short current_intensity = 0;
       boolean flag = false;
       TFColor basic_color = new TFColor(); 
      TFColor aug_color = new TFColor();
        aug_color.a = 0;
       aug_color.r = 0;
       aug_color.g = 0;
        aug_color.b = 0;

        VectorMath.setVector(entry_exit_vector, exitPoint[0]-entryPoint[0], exitPoint[1]-entryPoint[1], exitPoint[2]-entryPoint[2]);
        for (double current_dis = 0; current_dis < total_dis && aug_color.a<1; current_dis += sampleStep){
            for(int i = 0; i<3; i++){
                current_point[i] = (current_dis/total_dis)*entry_exit_vector[i]+entryPoint[i];
            }
            current_intensity = volume.getVoxelInterpolate(current_point);
            basic_color = tFunc.getColor(current_intensity);
            float alpha = (float) (1 - Math.pow(1-basic_color.a, sampleStep));
            aug_color.r += (1-aug_color.a) * basic_color.r *alpha;
            aug_color.g += (1-aug_color.a) * basic_color.g *alpha;
            aug_color.b += (1-aug_color.a) * basic_color.b *alpha;
            aug_color.a += (1-aug_color.a) * alpha;
            flag = (current_dis+sampleStep) >= total_dis ? true :false;
        }
        if (flag){
            current_intensity = volume.getVoxelInterpolate(exitPoint);
            basic_color = tFunc.getColor(current_intensity);
           float alpha = (float) (1 - Math.pow(1-basic_color.a, sampleStep));
           aug_color.r += (1-aug_color.a) * basic_color.r *alpha;
           aug_color.g += (1-aug_color.a) * basic_color.g *alpha;
           aug_color.b += (1-aug_color.a) * basic_color.b *alpha;
           aug_color.a += (1-aug_color.a) * alpha;
       }
       return doublesToColor(1,aug_color.r,aug_color.g,aug_color.b);
   }
    
   
   /*----------------------*/
      //2D Transfer fuction with shading 
   /*----------------------*/
    private void traceRay2Dtransferfucntion(double[] viewMatrix){

        // clear image        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // image is square
        int imageCenter = image.getWidth() >> 1;

        double[] voxelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() >> 1, volume.getDimY() >> 1, volume.getDimZ() >> 1);

        // sample on a plane through the origin of the volume data
        TFColor voxelColor = new TFColor();
        short voxelIntensity = 0;
        VoxelGradient voxelGradient = new VoxelGradient();
        
        // parameters to be used to shade
        double [] NormVec = new double[3];
        double dotProductNL = 0;
        double dotProductNH = 0;
        double colorIncrement = 0;
        
        
        // interactive mode
        int sampleStep = 1;
        if(interactiveMode == true) {
            sampleStep = 5;
        }
        
        // get base intensity, radius, and base color from TransferFunction2D
        final short baseIntensity = tfEditor2D.triangleWidget.baseIntensity;
        final double radius = tfEditor2D.triangleWidget.radius;
        final TFColor baseColor = tfEditor2D.triangleWidget.color;
        
        // extend widget
        final double lowGradientMagnitude = tfEditor2D.triangleWidget.minGradient;
        final double upGradientMagnitude = tfEditor2D.triangleWidget.maxGradient;

        final double XStep = viewVec[0] * sampleStep;
        final double YStep = viewVec[1] * sampleStep;
        final double ZStep = viewVec[2] * sampleStep;        
        
        for (int j = 0; j < image.getHeight(); j ++) {
            
            double voxelCoordXStart = uVec[0] * (-imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
            double voxelCoordYStart = uVec[1] * (-imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
            double voxelCoordZStart = uVec[2] * (-imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];
            
            for (int i = 0; i < image.getWidth(); i++) {
                
                TFColor pixelColor = new TFColor(0, 0, 0, 1);
                                
                // compute intersections
                long tXMin = Long.MIN_VALUE;
                long tXMax = Long.MAX_VALUE;
                long tYMin = Long.MIN_VALUE;
                long tYMax = Long.MAX_VALUE;
                long tZMin = Long.MIN_VALUE;
                long tZMax = Long.MAX_VALUE;
                boolean isIntersected = true;
                
                if(viewVec[0] != 0) {
                    tXMin = Math.round(-voxelCoordXStart / viewVec[0]);
                    tXMax = Math.round((volume.getDimX() - voxelCoordXStart) / viewVec[0]);;
                }
                else {
                    if(voxelCoordXStart < 0 || voxelCoordXStart >= volume.getDimX()) {
                        isIntersected = false;
                    }
                }
                
                if(viewVec[1] != 0) {
                    tYMin = Math.round(-voxelCoordYStart / viewVec[1]);
                    tYMax = Math.round((volume.getDimY() - voxelCoordYStart) / viewVec[1]);
                }
                else {
                    if(voxelCoordYStart < 0 || voxelCoordYStart >= volume.getDimY()) {
                        isIntersected = false;
                    }
                }
                
                if(viewVec[2] != 0) {
                    tZMin = Math.round(-voxelCoordZStart / viewVec[2]);
                    tZMax = Math.round((volume.getDimZ() - voxelCoordZStart) / viewVec[2]);
                }
                else {
                    if(voxelCoordZStart < 0 || voxelCoordZStart >= volume.getDimZ()) {
                        isIntersected = false;
                    }
                }
                
                //swap value if min larger than max
                if(tXMin > tXMax) {
                    tXMin = tXMin + tXMax;
                    tXMax = tXMin - tXMax;
                    tXMin = tXMin - tXMax;
                }
                //swap value if min larger than max
                if(tYMin > tYMax) {
                    tYMin = tYMin + tYMax;
                    tYMax = tYMin - tYMax;
                    tYMin = tYMin - tYMax;          
                }
                //swap value if min larger than max
                if(tZMin > tZMax) {
                    tZMin = tZMin + tZMax;
                    tZMax = tZMin - tZMax;
                    tZMin = tZMin - tZMax;          
                }
                
                long start = Math.max(tXMin, Math.max(tYMin, tZMin));
                long end = Math.min(tXMax, Math.min(tYMax, tZMax));

                if(start > end) {
                    isIntersected = false;
                }
                
                // ray intersects with volume
                if(isIntersected) {
                    voxelCoord[0] = voxelCoordXStart + (end + sampleStep) * viewVec[0];
                    voxelCoord[1] = voxelCoordYStart + (end + sampleStep) * viewVec[1];
                    voxelCoord[2] = voxelCoordZStart + (end + sampleStep) * viewVec[2];
                    
                    for(long u = end; u > start; u -= sampleStep) {
                        
                        // move to the next sample
                        voxelCoord[0] -= XStep;
                        voxelCoord[1] -= YStep;
                        voxelCoord[2] -= ZStep;

                        // get voxel intensity
                        voxelIntensity =volume.getVoxelInterpolate(voxelCoord);                        
                        voxelGradient = gradients.getGradient(voxelCoord);

                        voxelColor.r = baseColor.r;
                        voxelColor.g = baseColor.g;
                        voxelColor.b = baseColor.b;

                        // re-weigh the opacity
                        double absDiffGradRatio = Math.abs(voxelIntensity - baseIntensity) / (voxelGradient.mag + 1e-6);

                        if(voxelGradient.mag < lowGradientMagnitude || voxelGradient.mag > upGradientMagnitude) {
                            continue;
                        }
                        
                        // opacity reweighting
                        if(voxelGradient.mag <= 1e-6 && voxelIntensity == baseIntensity) {
                            voxelColor.a = baseColor.a;
                        }
                        else if(voxelGradient.mag > 1e-6 && absDiffGradRatio <= radius) {
                            voxelColor.a = baseColor.a * (1.0 - 1.0 / radius * absDiffGradRatio);
                        }
                        else{
                            continue;
                        }
                        
                        // shading
                        if(shadingMode) {
                            // surface normal at voxel
                            NormVec[0] = voxelGradient.x;
                            NormVec[1] = voxelGradient.y;
                            NormVec[2] = voxelGradient.z;
                            dotProductNL = Math.max(VectorMath.dotproduct(viewVec, NormVec) / (VectorMath.length(NormVec) + 1e-6), 0);
                            dotProductNH = dotProductNL;
                            colorIncrement = SHADING_SPEC_COEFF * Math.pow(dotProductNH, SHADING_ALPHA);
                            voxelColor.r = SHADING_AMBIENT_COEFF.r + voxelColor.r * (SHADING_DIFF_COEFF * dotProductNL) + colorIncrement;
                            voxelColor.g = SHADING_AMBIENT_COEFF.g + voxelColor.g * (SHADING_DIFF_COEFF * dotProductNL) + colorIncrement;
                            voxelColor.b = SHADING_AMBIENT_COEFF.b + voxelColor.b * (SHADING_DIFF_COEFF * dotProductNL) + colorIncrement;
                        }

                        // composite voxel colors
                        pixelColor.r = (1 - voxelColor.a) * pixelColor.r + voxelColor.a * voxelColor.r;
                        pixelColor.g = (1 - voxelColor.a) * pixelColor.g + voxelColor.a * voxelColor.g;
                        pixelColor.b = (1 - voxelColor.a) * pixelColor.b + voxelColor.a * voxelColor.b;
                        pixelColor.a = (1 - voxelColor.a) * pixelColor.a + voxelColor.a;
                    }
                                        
                    // BufferedImage expects a pixel color packed as ARGB in an int
                    int c_alpha = pixelColor.a <= 1.0 ? (int) Math.floor(pixelColor.a * 255) : 255;
                    int c_red = pixelColor.r <= 1.0 ? (int) Math.floor(pixelColor.r * 255) : 255;
                    int c_green = pixelColor.g <= 1.0 ? (int) Math.floor(pixelColor.g * 255) : 255;
                    int c_blue = pixelColor.b <= 1.0 ? (int) Math.floor(pixelColor.b * 255) : 255;
                    int finalPixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                    image.setRGB(i, j, finalPixelColor);
                }
                else {
                    // background
                    image.setRGB(i, j, 1 << 24);
                }
                // move to the next pixel
                voxelCoordXStart += uVec[0];
                voxelCoordYStart += uVec[1];
                voxelCoordZStart += uVec[2];
            } //end imageWidth
        } // end imageHeight
    }

    /*------------------------*/
    /*for converting color*/
        private int doublesToColor(double a, double r, double g, double b) {
        int c_alpha = a <= 1.0 ? (int) Math.floor(a * 255) : 255;
        int c_red = r <= 1.0 ? (int) Math.floor(r * 255) : 255;
        int c_green = g <= 1.0 ? (int) Math.floor(g * 255) : 255;
        int c_blue = b <= 1.0 ? (int) Math.floor(b * 255) : 255;
        int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
        return pixelColor;
    }
    /*---------------*/
            
     /*---------------*/
        //MIP
  int traceRayMIP(double[] entryPoint, double[] exitPoint, double[] viewVec, double sampleStep) {
        /* to be implemented:  You need to sample the ray and implement the MIP
         * right now it just returns yellow as a color
        */
        
        short max = 0; 
        int direction = direction(entryPoint, exitPoint);
        sampleStep *= direction;
         // interactive mode
     //   int sampleStep = 1;
        if(interactiveMode == true) {
            sampleStep = 5;
        }
        if (direction < 0){
            double[] position = entryPoint;
            while (between(position, entryPoint, exitPoint)) {
                max = (short)Math.max(max, volume.getVoxelInterpolate(position));
             //   System.out.println("samplestep is"+sampleStep);
                position = add(position, scale(sampleStep, viewVec));
            }
        }
        else{
            double[] position = exitPoint;
            while (between(position, exitPoint, entryPoint)) {
            max = (short)Math.max(max, volume.getVoxelInterpolate(position));
            position = add(position, scale(sampleStep, viewVec));
            } 
        }
         /*Applying TF color to MIP*/
      /*   TFColor auxColor = new TFColor(); 
        auxColor = tFunc.getColor(max);
        voxelColor.r=auxColor.r;voxelColor.g=auxColor.g;voxelColor.b=auxColor.b;voxelColor.a=auxColor.a;
        
         // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int color = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
      return color;*/
        return (max << 24) | (255 << 16) | (255 << 8) | 255;
    }
  boolean between(double[] position, double[] a, double[] b) {
        for (int i = 0; i < position.length; i += 1) {
            if (position[i]  > a[i] && position[i] > b[i] || position[i] < a[i] && position[i] < b[i]) {
                return false;
            }
        }
        return true;
    }
    
    /*
     * Computes whether a is "before" b (-1) or after b (+1).
    */
    int direction(double[] a, double[] b) {
        for (int i = 0; i < a.length; i += 1) {
            if (a[i] > b[i]) {
                return -1;
            }
        }
        return 1;
    }
    
     /**
     * Computes the multiplication of the a vector by the factor s.
     */
    double[] scale(double s, double[] a) {
        double[] b = new double[a.length];
        for (int i = 0; i < a.length; i += 1) {
            b[i] = s * a[i];
        }
        return b;
    }
    
    /*
     * Computes the element-by-element addition of vectors a and b.
    */
    double[] add(double[] a, double[] b) {
        double[] c = new double[a.length];
        for (int i = 0; i < a.length; i += 1) {
            c[i] = a[i] + b[i];
        }
        return c;
    }
    
    /*-----------------*/
    
   //Entry point and exit point computation
    void computeEntryAndExit(double[] p, double[] viewVec, double[] entryPoint, double[] exitPoint) {

        for (int i = 0; i < 3; i++) {
            entryPoint[i] = -1;
            exitPoint[i] = -1;
        }

        double[] plane_pos = new double[3];
        double[] plane_normal = new double[3];
        double[] intersection = new double[3];

        VectorMath.setVector(plane_pos, volume.getDimX(), 0, 0);
        VectorMath.setVector(plane_normal, 1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, -1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, volume.getDimY(), 0);
        VectorMath.setVector(plane_normal, 0, 1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, -1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, volume.getDimZ());
        VectorMath.setVector(plane_normal, 0, 0, 1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, 0, -1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

    }

    /*--------------------*/
    //For MIP, compositing and 2D transfer functions
    /*--------------------*/
    void raycast(double[] viewMatrix) {
        /* To be partially implemented:
            This function traces the rays through the volume. Have a look and check that you understand how it works.
            You need to introduce here the different modalities MIP/Compositing/TF2/ etc...*/
        
        //Pass view matrix directly for 2d transfer function (Shading calculation) 
        if(tf2dMode)
        {
            
            traceRay2Dtransferfucntion(viewMatrix);
        }
        
        //for mip and compositing
        else
        {
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);


        int imageCenter = image.getWidth() / 2;
        
       //double limit = 1.0;
        
               
        double[] pixelCoord = new double[3];
        double[] entryPoint = new double[3];
        double[] exitPoint = new double[3];
        int increment=1;
        float sampleStep=0.2f;
        

        System.out.println(viewVec[0] + " " + viewVec[1] + " " + viewVec[2]);
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        

        for (int j = 0; j < image.getHeight(); j += increment) {
            for (int i = 0; i < image.getWidth(); i += increment) {
                // compute starting points of rays in a plane shifted backwards to a position behind the data set
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) - viewVec[0] * imageCenter
                        + volume.getDimX() / 2.0;
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) - viewVec[1] * imageCenter
                        + volume.getDimY() / 2.0;
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) - viewVec[2] * imageCenter
                        + volume.getDimZ() / 2.0;

                computeEntryAndExit(pixelCoord, viewVec, entryPoint, exitPoint);
                if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                  
                   int pixelColor = 0;
                                   
                    /* set color to green if MipMode- see slicer function*/
                   if(mipMode && !super.interactiveMode) 
                        pixelColor= traceRayMIP(entryPoint,exitPoint,viewVec,sampleStep);
                   if(compositingMode && !super.interactiveMode)
                         pixelColor= traceRayCompositingF2B(entryPoint,exitPoint,viewVec,sampleStep);
                   
                    for (int ii = i; ii < i + increment; ii++) {
                        for (int jj = j; jj < j + increment; jj++) {
                            image.setRGB(ii, jj, pixelColor);
                        }
                    }
                }
             
            }
         }
      }

    }

    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0]; // how far is it from the origin multiplied by the perp vector. xdim
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1]; // ydim
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2]; //zdim

                int val = volume.getVoxelInterpolate(pixelCoord);
                // Map the intensity to a grey value by linear scaling
               /* voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
               voxelColor.a = val > 0 ? 1.0 : 0.0;  */ // this makes intensity 0 completely transparent and the rest opaque
                
                // Alternatively, apply the transfer function to obtain a color
                TFColor auxColor = new TFColor(); 
                auxColor = tFunc.getColor(val);
                voxelColor.r=auxColor.r;voxelColor.g=auxColor.g;voxelColor.b=auxColor.b;voxelColor.a=auxColor.a;
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
    }


    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        if (slicerMode) {
            slicer(viewMatrix);    
        } else {
            raycast(viewMatrix);
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}

