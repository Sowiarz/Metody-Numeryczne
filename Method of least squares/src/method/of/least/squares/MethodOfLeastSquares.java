/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package method.of.least.squares;
import Jama.*;
import java.awt.geom.Point2D;
/**
 *
 * @author Sowiarz
 */
public class MethodOfLeastSquares {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        double[][] p = {{1.,1.},{2.,3.},{3.,4.},{4.,5.},{5.,3},{6.,4.}};
        Point2D[] points = {
                new Point2D.Double(1,1),
                new Point2D.Double(2,3),
                new Point2D.Double(3,4),
                new Point2D.Double(4,5),
                new Point2D.Double(5,3),
                new Point2D.Double(6,4)
        };
        Matrix point = new Matrix(p);
        int degree = 1;
        
        Matrix A = createA(points, degree);
        Matrix B = createB(points);

        Matrix AT = A.transpose();
        Matrix ATimesAT = AT.times(A);
        ATimesAT.print(6, 3);
        Matrix ATTimesB = AT.times(B);
        ATTimesB.print(2, 2);
        Matrix test = ATimesAT.solve(ATTimesB);
        test.print(6, 3);
       
       
       int m = point.getRowDimension();
       int n = point.getColumnDimension();
       
       Matrix Q = new Matrix(m,n); 
       Matrix R = new Matrix(n,n); 
            
       for (int k=0 ; k < n ; k++) {
         Q.setMatrix(0, m-1, k, k, point.getMatrix(0, m-1, k, k));
         for (int i=0 ; i < k ; i++) {
             double rik=(Q.getMatrix(0, m-1, i, i).transpose().  
             times(Q.getMatrix(0, m-1, k, k))).get(0,0);  
             R.set(i, k, rik);
             Q.setMatrix(0, m-1,k,k,Q.getMatrix(0, m-1, k, k).
                     minus(Q.getMatrix(0, m-1, i, i).times(rik)));
         }
         double rkk =Q.getMatrix(0, m-1, k, k).norm2();
         R.set(k, k, rkk);
         Q.setMatrix(0,m-1,k,k,Q.getMatrix(0,m-1,k,k).times(1/rkk));
       }
       System.out.println("Macierz Q :");
       Q.print(10,6);
       System.out.println("Macierz R :");
       R.print(10,6);
       Q.times(R).print(10,6);
       Matrix Qt = Q.transpose();
       Matrix QtB = Qt.times(B);
       R.solve(QtB);
       System.out.printf("WYNIK\n");
       R.print(10, 6);
       int i,j;
       System.out.printf("Aproksymacja wielomianowa BIBLIOTECZNA\n");
       
       int rowc=point.getRowDimension();
       System.out.printf("Zbiór punktów do aproksymacji : \n");
       point.print(5,5);
       Matrix Y=point.getMatrix(0,rowc-1,1,1);
       Matrix X=new Matrix(rowc,degree);
       
       for ( i=0 ; i< rowc ; i++){
        for( j=0 ; j< degree ; j++){
           if(j==0) {X.set(i, j, 1.0);} else
                    {X.set(i, j,Math.pow(A.get(i, 0),j));}  
        }  
       }
       System.out.println("Macierz pomocnicza :");
       X.print(10,3);
       //Współczynniki wielomianu :
       QRDecomposition QR=X.qr();
       Matrix Qb = QR.getQ();
       Matrix Rb = QR.getR(); 
       Matrix wynik=R.solve(Q.transpose().times(Y));
       System.out.printf("Współczynniki wielomianu stopnia %d : \n",degree);
       wynik.print(1, 15);
    }
    public static Matrix createA(Point2D[] x, int degree){
        Matrix result = new Matrix(x.length, degree + 1);

        for (int i = 0; i < result.getRowDimension(); i++) {

            for (int j = 0; j < degree + 1; j++) {

                if(j == 0) result.set(i, 0, 1);
                else result.set(i, j, Math.pow(x[i].getX(), j));

            }
        }
        return result;
    }
    public static Matrix createB(Point2D[] x) {

        Matrix result = new Matrix(x.length, 1);

        for (int i = 0; i < result.getRowDimension(); i++) {

            result.set(i, 0, x[i].getY());

        }
        return result;

    }
    public static double horner(double[] factorTable, int degree, double x){
        if(degree==0)
            return factorTable[0];
        return x*horner(factorTable, degree-1, x)+factorTable[degree];
    }
    
}
