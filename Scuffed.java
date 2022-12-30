import java.util.*;
import java.io.*;

public class Scuffed {
    static int n, m;
    static int startI, startJ; // Stores starting position of each floodfill
    static int[][] course; // Stores the skiing course heights
    static boolean[][] waypoints; // Stored as booleans instead of 1s and 0s
    static boolean[][] vis; // Visited array for floodfill

    public static void main(String[] args) throws IOException {
        Kattio io = new Kattio("ccski");

        n = io.nextInt();
        m = io.nextInt();

        int minHeight = Integer.MAX_VALUE;
        int maxHeight = Integer.MIN_VALUE;

        course = new int[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                course[i][j] = io.nextInt();
                minHeight = Math.min(minHeight, course[i][j]);
                maxHeight = Math.max(maxHeight, course[i][j]);
            }
        }

        waypoints = new boolean[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (io.nextInt() == 1) {
                    waypoints[i][j] = true;

                    // Keep one of the waypoints as the starting position
                    startI = i;
                    startJ = j;
                } else {
                    waypoints[i][j] = false;
                }
            }
        }

        // Binary search for the minimum possible value of d that works
        // We can set "hi" to maxHeight - minHeight because that will be the
        // maximum value of d we need to try
        int lo = 0;
        int hi = maxHeight - minHeight;
        int minD = -1;
        while (lo <= hi) {
            int d = (lo + hi) / 2;
            if (reachable(d)) {
                minD = d;
                hi = d - 1;
            } else {
                lo = d + 1;
            }
        }

        io.println(minD);
        io.close();
    }

    // i and j store current position
    // d stores current value of d
    // prevHeight stores the height of the previous cell
    static void floodfill(int i, int j, int d, int prevHeight) {
        // Check if we are out of bounds
        if (i < 0 || i >= n || j < 0 || j >= m) {
            return;
        }
        // Check if current position has already been visited
        if (vis[i][j]) {
            return;
        }
        // Check if current position can be visited from previous position
        if (Math.abs(course[i][j] - prevHeight) > d) {
            return;
        }

        // Mark position as visited if all checks pass
        vis[i][j] = true;

        // Visit each adjacent cell
        floodfill(i + 1, j, d, course[i][j]);
        floodfill(i - 1, j, d, course[i][j]);
        floodfill(i, j + 1, d, course[i][j]);
        floodfill(i, j - 1, d, course[i][j]);
    }

    static boolean reachable(int d) {
        // Reset visited array and begin floodfill (DFS) from start position
        vis = new boolean[n][m];
        floodfill(startI, startJ, d, course[startI][startJ]);

        // Check each cell: if it is a waypoint and it hasn't been visited,
        // we know not all waypoints are reachable from one another
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (waypoints[i][j] && !vis[i][j]) {
                    return false;
                }
            }
        }

        return true;
    }

    static class Kattio extends PrintWriter {
        private BufferedReader r;
        private StringTokenizer st;

        // standard input
        public Kattio() {
            this(System.in, System.out);
        }

        public Kattio(InputStream i, OutputStream o) {
            super(o);
            r = new BufferedReader(new InputStreamReader(i));
        }

        // USACO-style file input
        public Kattio(String problemName) throws IOException {
            super(new FileWriter(problemName + ".out"));
            r = new BufferedReader(new FileReader(problemName + ".in"));
        }

        // returns null if no more input
        public String next() {
            try {
                while (st == null || !st.hasMoreTokens())
                    st = new StringTokenizer(r.readLine());
                return st.nextToken();
            } catch (Exception e) {
            }
            return null;
        }

        public int nextInt() {
            return Integer.parseInt(next());
        }

        public double nextDouble() {
            return Double.parseDouble(next());
        }

        public long nextLong() {
            return Long.parseLong(next());
        }

        public String readLine() throws Exception {
            return r.readLine();
        }

        public long[] readArrayLong(int n) {
            long[] arr = new long[n];
            for (int i = 0; i < n; i++)
                arr[i] = nextLong();
            return arr;
        }

        public Integer[] readArrayI(int n) {
            Integer[] arr = new Integer[n];
            for (int i = 0; i < n; i++)
                arr[i] = nextInt();
            return arr;
        }

        public int[] readArray(int n) {
            int[] arr = new int[n];
            for (int i = 0; i < n; i++)
                arr[i] = nextInt();
            return arr;
        }

        public int[]prefixSum(int N, int[]a){
            int[]nw = new int[N + 1];
            nw[0] = 0;
            for (int i = 0; i < N; i ++) {
                nw[i + 1] = nw[i] + a[i];
            }

            return nw;
        }

        public int[][] prefixSum(int N, int[][]grid){
            int[][]pfx = new int[grid.length][grid.length];
            for (int i = 1; i < N + 1; i++) {
                for (int j = 1; j < N + 1; j++) {
                    pfx[i][j] = grid[i][j]
                            + pfx[i-1][j]
                            + pfx[i][j-1]
                            - pfx[i-1][j-1];
                }
            }

            return pfx;
        }
        public int prefixSum2D(int x1, int y1, int x2, int y2, int[][]pfx) {
            return pfx[x2][y2] - pfx[x1-1][y2] - pfx[x2][y1-1] + pfx[x1-1][y1-1];
        }
    }
}