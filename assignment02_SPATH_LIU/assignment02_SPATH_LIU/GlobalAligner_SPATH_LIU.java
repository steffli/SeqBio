package assignment02;

import assignment01.FastA_SPATH_LIU;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 * GlobalAligner_SPATH_LIU
 * Sequence Bioinformatics, WS 22/23
 */
public class GlobalAligner_SPATH_LIU {
	// variables for the algorithms
	static int match = 1;
	static int mismatch = -1;
	static int d = 1;
	public static ArrayList<ArrayList<Integer>> tracebackLin = new ArrayList<ArrayList<Integer>>();


	public static void main(String[] args) throws IOException {
		if(args.length<1 || args.length>2)
			throw new IOException("Usage: GlobalAligner_SPATH_LIU infile [quadraticSpace|linearSpace|noDP]");

		var list=FastA_SPATH_LIU.read(args[0]);

		if(list.size()!=2)
			throw new IOException("Wrong number of input sequences: "+list.size());

		var mode=(args.length==2?args[1]:"quadraticSpace");


		switch(mode) {
			case "quadraticSpace" -> runNeedlemanWunschQuadraticSpace(list.get(0),list.get(1));
			case "linearSpace" ->runNeedlemanWunschLinearSpace(list.get(0),list.get(1));
			case "noDP" ->runNeedlemanWunschRecursively(list.get(0),list.get(1));
			default -> throw new IOException("Unknown mode: "+mode);
		}

	}

	/**
	 * computes the optimal global alignment score and an alignment, using quadratic space.
	 * Prints out the optimal score and a corresponding alignment
	 * Also prints out the number of milliseconds required for the computation
	 * @param x Pair for computation
	 * @param y Pair for computation
	 */
	public static void runNeedlemanWunschQuadraticSpace(FastA_SPATH_LIU.Pair x, FastA_SPATH_LIU.Pair y) {
		// todo: implement, Assignment 2.1
		String s1 = x.sequence();
		String s2 = y.sequence();

		// initialization
		Integer[][] matrix = new Integer[s1.length()+1][s2.length()+1];
		Integer[][] traceback = new Integer[s1.length()+1][s2.length()+1];
		// 0 means up, 1 left, 2 match/mismatch, 3 stop


		long start = System.nanoTime();


		traceback[0][0] = 3;
		for (int i = 1; i <= s1.length(); i++) {
			traceback[i][0] = 1;
		}
		for (int j = 1; j <= s2.length(); j++) {
			traceback[0][j] = 0;
		}
		for (int i = 0; i <= s1.length(); i++) {
			matrix[i][0] = -i * d;
		}
		for (int j = 0; j <= s2.length(); j++) {
			matrix[0][j] = -j * d;
		}
		// recursion
		int score;
		for (int i = 1; i < matrix.length; i++) {
			for (int j = 1; j < matrix[i].length; j++) {
				if (s1.charAt(i-1) == s2.charAt(j-1)) score = match;
				else score = mismatch;

				int m = matrix[i-1][j-1]+score;
				int l = matrix[i-1][j]-d;
				int u = matrix[i][j-1]-d;
				if ((m>l) && (m>u)) {
					matrix[i][j] = m;
					traceback[i][j] = 2;
				}
				else if (l>u) {
					matrix[i][j] = l;
					traceback[i][j] = 1;
				}
				else {
					matrix[i][j] = u;
					traceback[i][j] = 0;
				}
				// matrix[i][j] = Math.max(matrix[i-1][j-1]+score, Math.max(matrix[i-1][j]-d, matrix[i][j-1]-d));

			}
		}
		System.out.println("optimal score is: " + matrix[matrix.length - 1][matrix[0].length - 1]);
		// traceback
		StringBuilder align1 = new StringBuilder();
		StringBuilder align2 = new StringBuilder();

		int i = matrix.length - 1;
		int j = matrix[0].length - 1;
		int tb;
		do {
			tb = traceback[i][j];
			if (tb == 2) {
				i -= 1;
				j -= 1;
				align1.append(s1.charAt(i));
				align2.append(s2.charAt(j));
			} else if (tb == 1) {
				i -= 1;
				align1.append(s1.charAt(i));
				align2.append('-');
			} else if (tb == 0) {
				j -= 1;
				align2.append(s2.charAt(j));
				align1.append('-');
			}

		} while (tb != 3);

		// measure time and transform to milliseconds
		long end = System.nanoTime();

		System.out.println(x.header() + ":\t" + align1.reverse());
		System.out.println(y.header() + ":\t" + align2.reverse());
		System.out.println("milliseconds used for computation: " + (end - start) / 1000000);


	}

	/**
	 * computes the optimal global alignment score and an alignment, using linear space
	 * Prints out the optimal score and a corresponding alignment
	 * Also prints out the number of milliseconds required for the computation
	 * @param x
	 * @param y
	 */
	public static void runNeedlemanWunschLinearSpace(FastA_SPATH_LIU.Pair x, FastA_SPATH_LIU.Pair y) {
		// todo: implement, Assignment 2.2
		String s1 = x.sequence();
		String s2 = y.sequence();
		// initialization
		Integer[] array1 = new Integer[s2.length() + 1];
		Integer[] array2 = new Integer[s2.length() + 1];
		Integer[] tb1 = new Integer[s2.length() + 1];
		Integer[] tb2 = new Integer[s2.length() + 1];

		long start = System.nanoTime();


		for (int j = 0; j <= s2.length(); j++) {
			array1[j] = -j * d;
			tb1[j] = j;
		}
		// recursion
		int score;
		for (int i = 1; i < s1.length() + 1; i++) {
			array2[0] = -i * d;
			for (int j = 1; j < array1.length; j++) {
				if (s1.charAt(i-1) == s2.charAt(j-1)) score = match;
				else score = mismatch;

				int m = array1[j-1]+score;
				int l = array1[j]-d;
				int u = array2[j-1]-d;
				if ((m>l) && (m>u)) {
					array2[j] = m;
					if(i>((s1.length()+1)/2)) tb2[j] = tb1[j-1];
				}
				else if (l>u) {
					array2[j] = l;
					if(i>((s1.length()+1)/2)) tb2[j] = tb1[j];
				}
				else {
					array2[j] = u;
					if(i>((s1.length()+1)/2)) tb2[j] = tb2[j-1];
				}
				// array2[j] = Math.max(array1[j-1]+score, Math.max(array1[j]-d, array2[j-1]-d));
			}
			Integer[] swap = array1;
			array1 = array2;
			array2 = swap;
			if(i>((s1.length()+1)/2)) {
				Integer[] swaptb = tb1;
				tb1 = tb2;
				tb2 = swaptb;
			}
		}
		// array1 always last updated
		int optimalScore = array1[array1.length-1];
		int c =  tb1[tb1.length-1];
		/*if (s1.length() % 2 != 0) {
			optimalScore = array1[array1.length-1];
			c = tb1[tb1.length-1];
		}
		else {
			optimalScore = array2[array2.length-1];
			c = tb2[tb2.length-1];
		}*/
		System.out.println("optimal score is: " + optimalScore);
		System.out.println("first c is column "+(s1.length() + 1)/2+ " row " + c);

		// traceback

		// save rows of traceback to arraylist(columns are indices)
		for (int i = 0; i<s2.length() + 1; i++) {
			tracebackLin.add(new ArrayList<>());
		}

		int midRow =  tb1[tb1.length-1];
		int midCol = (s1.length()+1)/2;

		// add mid field to traceback
		ArrayList<Integer> row = tracebackLin.get(midRow);
		row.add(midCol);
		tracebackLin.set(midRow, row);
		//add last field to traceback
		row = tracebackLin.get(s2.length());
		row.add(s1.length());
		tracebackLin.set(s2.length(), row);
		//add first field to traceback
		/*row = tracebackLin.get(0);
		row.add(0);
		tracebackLin.set(0, row);*/

		// compute c recursively
		computeC(0,0, s1.substring(0,midCol), s2.substring(0,midRow));
		computeC(midRow, midCol, s1.substring(midCol), s2.substring(midRow));

		// wie oben computen, und jeweils neu initialisieren, vielleicht extra Methode?


		// measure time and transform to milliseconds
		long end = System.nanoTime();
		System.out.println("milliseconds used for computation: " + (end - start) / 1000000);

		// write to console

		StringBuilder align1 = new StringBuilder();
		StringBuilder align2 = new StringBuilder();

		int prevCol = 0;
		int prevRow = 0;

		for (int a = 0; a < tracebackLin.size(); a++) {
			Collections.sort(tracebackLin.get(a));
			for (int b : tracebackLin.get(a)) {
				if (b>prevCol && a == prevRow) {
					align1.append(s1.charAt(b-1));
					align2.append("-");
				}
				else if (b==prevCol && a > prevRow) {
					align1.append("-");
					align2.append(s2.charAt(a-1));
				}
				else {
					align1.append(s1.charAt(b-1));
					align2.append(s2.charAt(a-1));
				}
				prevCol = b;
				prevRow = a;
			}
		}

		System.out.println(align1);
		System.out.println(align2);

		for (int a = 0; a < tracebackLin.size(); a++) {
			int count = 0;
			System.out.println("row " + a + " columns: ");
			for (int b : tracebackLin.get(a)) {
				System.out.print(b + " ");
			}
			System.out.println();
		}


	}

	public static void computeC(int realRow, int realColumn, String s1, String s2) {

		if (s1.length()<2 && s2.length()<2) return;
		if (s1.length()<2 || s2.length()<2) {
			if (s1.length() > s2.length()) {
				for (int i = 0; i < s1.length()-1; i++) {
					ArrayList<Integer> row = tracebackLin.get(realRow);
					row.add(realColumn+i);
					tracebackLin.set(realRow, row);
				}
			}
			else {
				for (int i = 1; i < s2.length(); i++) {
					ArrayList<Integer> row = tracebackLin.get(realRow+i);
					row.add(realColumn);
					tracebackLin.set(i+realRow, row);
				}
			}

			return;
		}

		// initialize arrays
		Integer[] array1 = new Integer[s2.length() + 1];
		Integer[] array2 = new Integer[s2.length() + 1];
		Integer[] tb1 = new Integer[s2.length() + 1];
		Integer[] tb2 = new Integer[s2.length() + 1];

		for (int j = 0; j <= s2.length(); j++) {
			array1[j] = -j * d;
			tb1[j] = j;
		}

		int midCol;
		midCol = (s1.length()+1)/2;

		int score;
		for (int i = 1; i < s1.length() + 1; i++) {
			array2[0] = -i * d;
			for (int j = 1; j < array1.length; j++) {
				if (s1.charAt(i-1) == s2.charAt(j-1)) score = match;
				else score = mismatch;

				int m = array1[j-1]+score;
				int l = array1[j]-d;
				int u = array2[j-1]-d;
				if ((m>l) && (m>u)) {
					array2[j] = m;
					if(i>midCol) tb2[j] = tb1[j-1];
				}
				else if (l>u) {
					array2[j] = l;
					if(i>midCol) tb2[j] = tb1[j];
				}
				else {
					array2[j] = u;
					if(i>midCol) tb2[j] = tb2[j-1];
				}
			}
			Integer[] swap = array1;
			array1 = array2;
			array2 = swap;
			if(i>midCol) {
				Integer[] swaptb = tb1;
				tb1 = tb2;
				tb2 = swaptb;
			}
		}


		int midRow = 0;
		if (tb1[tb1.length-1] != null) midRow =  tb1[tb1.length-1];
		int realMidRow = midRow + realRow;
		int realMidCol = midCol + realColumn;

		ArrayList<Integer> row = tracebackLin.get(realMidRow);
		row.add(realMidCol);
		tracebackLin.set(realMidRow, row);


		computeC(realRow, realColumn, s1.substring(0,midCol), s2.substring(0,midRow));
		computeC(realMidRow, realMidCol, s1.substring(midCol), s2.substring(midRow));

	}

	/**
	 * computes the optimal global alignment score using a recursion and no table
	 * Prints out the optimal score
	 * Also prints out the number of milliseconds required for the computation
	 * @param x
	 * @param y
	 */
	public static void runNeedlemanWunschRecursively(FastA_SPATH_LIU.Pair x, FastA_SPATH_LIU.Pair y) {
		// todo: implement using recursive function computeF, , Assignment 2.3
		String s1 = x.sequence();
		String s2 = y.sequence();
		long start = System.nanoTime();

		System.out.println("optimal score is: " + computeF(s1.length(),s2.length(),s1,s2));
		// measure time and transform to milliseconds
		long end = System.nanoTime();
		System.out.println("milliseconds used for computation: " + (end - start) / 1000000);

	}

	public static int computeF(int i,int j, String s1, String s2) {
		// todo: implement
		if (i==0 && j==0) return 0;
		else if (i==0) return -j * d;
		else if (j==0) return -i * d;
		else {
			int score;
			if (s1.charAt(i-1) == s2.charAt(j-1)) score = match;
			else score = mismatch;
			return Math.max(computeF(i-1,j-1, s1, s2) + score,
					Math.max(computeF(i-1,j,s1,s2)-d, computeF(i,j-1,s1,s2)-d));
		}
	}
}
