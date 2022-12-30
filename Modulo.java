import java.io.*;
import java.util.*;

public class Modulo
{
    static InputReader in = new InputReader();
    static PrintWriter out = new PrintWriter(System.out);

    public static final int MN = 200020;

    static int N, M, ans;
    static int[] hd = new int[MN], nx = new int[MN], to = new int[MN], s = new int[MN], p = new int[MN];

    public static void adde(int u, int v, int id)
    {
        nx[id] = hd[u];
        hd[u] = id;
        to[id] = v;
    }
    public static void dfs(int node)
    {
        s[node]=1;
        for(int id=hd[node];id!=0;id=nx[id])
        {
            dfs(to[id]);
            s[node]+=s[to[id]];
        }
    }
    public static void main(String... args)
    {
        N=in.nextInt();
        for(int i=2;i<=N;++i)
        {
            p[i]=in.nextInt();
            adde(p[i], i, i);
        }
        dfs(1);
        for(int i=1;i<=N;++i)
        {
            out.print(s[i]-1);
            if(i<N) out.print(" ");
            else out.println();
        }
        out.close();
    }
    static class InputReader extends PrintWriter {
        private BufferedReader r;
        private StringTokenizer st;

        // standard input
        public InputReader() {
            this(System.in, System.out);
        }

        public InputReader(InputStream i, OutputStream o) {
            super(o);
            r = new BufferedReader(new InputStreamReader(i));
        }

        // USACO-style file input
        public InputReader(String problemName) throws IOException {
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