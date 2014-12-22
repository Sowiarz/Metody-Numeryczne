/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cholesky.decomposition;
import Jama.*;
/**
 *
 * @author Sowiarz
 */
public class CholeskyDecomposition {
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Matrix result;
        double[][] temp = {{2.,4.,6.},{4.,10.,20.},{6.,20.,52.}};
        Matrix example = new Matrix(temp);
        
        if(isSymetricMatrix(example)==true){
            System.out.println("Macierz jest symetryczna");
        }
        else System.out.println("Macierz nie jest symetrzyczna");
        if(isPositiveDefinite(example)==true){
            System.out.println("Macierz jest pozytywnie określona");
        }
        else System.out.println("Macierz nie jest pozytywnie określona");
        result = CholeskyDecomposition(example);
        showMatrix(result);
    }
    public static void showMatrix(Matrix x) {
        int n = x.getColumnDimension();
        
        for(int i=0; i<n; i++) {
            for(int j=0; j<n; j++) {
                System.out.print(x.get(i, j));
                System.out.print(" ");
            }
            System.out.println("");
        }
    }
    public static Matrix CholeskyDecomposition (Matrix x) {
 
     
      double[][] A = x.getArray();
      int n = x.getColumnDimension();
      double[][] L = new double[n][n];
      
      for (int j = 0; j < n; j++) {
         double[] Lrowj = L[j];
         double d = 0.0;
         for (int k = 0; k < j; k++) {
            double[] Lrowk = L[k];
            double s = 0.0;
            for (int i = 0; i < k; i++) {
               s += Lrowk[i]*Lrowj[i];
            }
            Lrowj[k] = s = (A[j][k] - s)/L[k][k];
            d = d + s*s;
         }
         d = A[j][j] - d;
         L[j][j] = Math.sqrt(Math.max(d,0.0));
         for (int k = j+1; k < n; k++) {
            L[j][k] = 0.0;
         }
      }
      return new Matrix(L,n,n);
   }
    public static boolean isPositiveDefinite(Matrix x){
        for(int i = 1; i<x.getRowDimension(); i++){
                double temp = x.getMatrix(0, i, 0, i).det();
                if(temp<=0){
                    return false;
                }
        }
        return true;
    }
    public static boolean isSymetricMatrix(Matrix x){
        double[][] A = x.getArray();
        int n = x.getColumnDimension();
        boolean check = checkDimension(x);
        if(check == true){
            for(int i = 0; i<n; i++){
                for(int j = 0; j<n; j++){
                    if(i==j) continue;
                    else
                        if(A[i][j]==A[j][i]) continue;
                        else return false;
                }
            }
            return true;
        }
        else 
            return false;        
    }
    public static boolean checkDimension(Matrix x){
        if(x.getColumnDimension()==x.getRowDimension())
            return true;
        else 
            return false;
    }
}
