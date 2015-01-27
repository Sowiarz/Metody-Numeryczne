/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lu.decomposition;
import Jama.*;

/**
 *
 * @author Sowiarz
 */
public class LUDecomposition {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Matrix A;
        //Matrix U = new Matrix();
        //Matrix L;
        double[][] temp = {{1.,5.,7.,6.},{-2.,3.,4.,8.},{13.,6.,2.,-3.},{3.,2.,6.,5.}};
        //double[][] tempa = {{1.,1.,1.,1.},{1.,1.,1.,1.},{1.,1.,1.,1.},{1.,1.,1.,1.}};
        A = new Matrix(temp);
        Matrix U = new Matrix(temp);
        Matrix L = Matrix.identity(A.getColumnDimension(), A.getRowDimension());
        
        for(int j = 1; j<A.getColumnDimension()-1; j++){
            L = getMatrixMiInv(L, j).times(L);
        }
        
       for(int i = 0; i<A.getColumnDimension()-1; i++){
            U = getMatrixMi(U, i).times(U);
        }
        
        U.print(5, 10);
        L.print(5, 10);
        L.times(U).print(5, 10);
        }
    public static Matrix getMatrixMi(Matrix A, int x){
        
        Matrix temp = new Matrix(A.getColumnDimension(), A.getRowDimension());
        for(int i = 0; i<A.getColumnDimension(); i++){
            for(int j = 0; j<A.getColumnDimension(); j++){
                if(i == j){
                    temp.set(i, i ,1);
                }
                if(i<j){
                    temp.set(i, j, 0);
                }
                if(j<i && j!=x){
                    temp.set(i, j, 0);
                }
                if(j<i && j == x){
                    double tmp = -(A.get(i, x))/(A.get(x, x));
                    temp.set(i, x, tmp);
                }
            }
        }
        return temp;
    }
    public static Matrix getMatrixMiInv(Matrix A, int x) {
        
        Matrix temp = new Matrix(A.getColumnDimension(), A.getRowDimension());
        for(int i = 0; i<A.getColumnDimension(); i++){
            for(int j = 0; j<A.getColumnDimension(); j++){
                if(i == j){
                    temp.set(i, i ,1);
                }
                if(i<j){
                    temp.set(i, j, 0);
                }
                if(j<i && j!=x){
                    temp.set(i, j, 0);
                }
                if(j<i && j == x){
                    double tmp = (A.get(i, x))/(A.get(x, x));
                    temp.set(i, x, tmp);
                }
            }
        }
        return temp;
    }
    
}
