import java.io.*;
import java.util.*;

public class Automate {
    public static void main(String[] args) throws Exception {
        //fastq_dump();
        //trimResults();
        //trimPairedResults();
        //alignReads();
        //alignPairedReads();
        //peakCall();
        //peakCallPairedEnd();
        //convertToCorrectFormat();
        //check();
        prune();
    }

    static void fastq_dump() throws Exception{
        ArrayList<String>testCases = new ArrayList<String>();
        /*
        testCases.add("SRR6493909"); testCases.add("SRR6493910"); testCases.add("SRR5280434"); testCases.add("SRR5280435");
        testCases.add("SRR7968891"); testCases.add("SRR7968892"); testCases.add("SRR7192349");
        testCases.add("SRR10820724"); testCases.add("SRR10820725");

         */
        testCases.add("SRR6233822"); testCases.add("SRR6233823"); testCases.add("SRR6233824"); testCases.add("SRR6233825");
        testCases.add("SRR6233826"); testCases.add("SRR6233827"); testCases.add("SRR6233828"); testCases.add("SRR6233829");

        for(String commands: testCases){
            String command = "fastq-dump --outdir " + commands + " --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ";
            command += "/Users/yuxiangzhang/Desktop/Pc3/" + commands +"/" + commands + ".sra";
            System.out.println(command);
            /*
            Process proc = Runtime.getRuntime().exec(command);


            BufferedReader reader = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            String line = "";
            while((line = reader.readLine()) != null) {
                System.out.print(line + "\n");
            }
            proc.waitFor();

             */
        }
    }
    static void trimResults() throws Exception {
        //single end reads
        ArrayList<String>testCases = new ArrayList<String>();
        /*
        testCases.add("SRR6493909"); testCases.add("SRR6493910"); testCases.add("SRR5280434"); testCases.add("SRR5280435");
        testCases.add("SRR7968891"); testCases.add("SRR7968892"); testCases.add("SRR7192349");
        testCases.add("SRR10820724"); testCases.add("SRR10820725");
         */

        testCases.add("SRR6233822"); testCases.add("SRR6233823"); testCases.add("SRR6233824"); testCases.add("SRR6233825");
        testCases.add("SRR6233826"); testCases.add("SRR6233827"); testCases.add("SRR6233828"); testCases.add("SRR6233829");
        for (String commands: testCases){
            String command = "java -jar trimmomatic-0.39.jar SE -phred33 /Users/yuxiangzhang/Desktop/Pc3/" + commands + "/" + commands + "_pass.fastq.gz " + commands + ".output.fastq.gz  ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
            System.out.println(command);
        }
    }
    static void trimPairedResults(){
        ArrayList<String>testCases = new ArrayList<String>();
        testCases.add("SRR5018853"); testCases.add("SRR7968966"); testCases.add("SRR7968967");
        String first = "java -jar trimmomatic-0.39.jar PE";
        String command = "/Users/yuxiangzhang/Desktop/\"Research paper\"/\"Paired End\"/";
        String end = "ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
        for (String test: testCases){
            String path = command + (test + "/" + test + "_pass_1.fastq.gz");
            String path2 = command + (test + "/" + test + "_pass_2.fastq.gz");
            //output
            String output1 = test + "p1.fastq.gz"; String output2 = test + "un1.fastq.gz";
            String output3 = test + "p2.fastq.gz"; String output4 = test + "un2.fastq.gz";

            String fullCommand = first + " " + path + " " + path2 + " " + output1 + " " + output2 + " " + output3 + " " + output4 + " " + end;
            System.out.println(fullCommand);
        }
    }
    static void alignReads()throws Exception {
        ArrayList<String>testCases = new ArrayList<String>();
        /*
        testCases.add("SRR6493909"); testCases.add("SRR6493910"); testCases.add("SRR5280434"); testCases.add("SRR5280435");
        testCases.add("SRR7968891"); testCases.add("SRR7968892"); testCases.add("SRR7192349");
        testCases.add("SRR10820724"); testCases.add("SRR10820725");
         */

        testCases.add("SRR6233822"); testCases.add("SRR6233823"); testCases.add("SRR6233824"); testCases.add("SRR6233825");
        testCases.add("SRR6233826"); testCases.add("SRR6233827"); testCases.add("SRR6233828"); testCases.add("SRR6233829");
        String commands = "bowtie2 --very-fast-local -x /Users/yuxiangzhang/Desktop/bowtie2_index/bowtie  -U /Users/yuxiangzhang/Desktop/Pc3/";
        for (String command: testCases){
            System.out.println(commands + command + "/" + command + ".output.fastq.gz -S /Users/yuxiangzhang/Desktop/Pc3/" + command + "/" + command + ".sam");
        }
    }
    static void alignPairedReads(){
        ArrayList<String>testCases = new ArrayList<String>();
        testCases.add("SRR5018853"); testCases.add("SRR7968966"); testCases.add("SRR7968967");
        String path = "bowtie2 --very-fast-local -x /Users/yuxiangzhang/Desktop/bowtie2_index/bowtie ";
        String prefix = " -S /Users/yuxiangzhang/Desktop/\"Research paper\"/\"Paired End\"/";
        String prefix1 = "-1 /Users/yuxiangzhang/Desktop/\"Research paper\"/\"Paired End\"/";
        String prefix2 = "-2 /Users/yuxiangzhang/Desktop/\"Research paper\"/\"Paired End\"/";
        for (String commands: testCases){
            System.out.println(path + prefix1 + commands + "/" + commands + "p1.fastq.gz " + prefix2 + commands + "/" + commands + "p2.fastq.gz" + prefix + commands + ".sam");
        }


    }
    static void peakCall(){
        ArrayList<String>testCases = new ArrayList<String>();
        /*
        testCases.add("SRR6493909"); testCases.add("SRR6493910"); testCases.add("SRR5280434"); testCases.add("SRR5280435");
        testCases.add("SRR7968891"); testCases.add("SRR7968892"); testCases.add("SRR7192349");
        testCases.add("SRR10820724"); testCases.add("SRR10820725");
        */
        testCases.add("SRR6233822"); testCases.add("SRR6233823"); testCases.add("SRR6233824"); testCases.add("SRR6233825");
        testCases.add("SRR6233826"); testCases.add("SRR6233827"); testCases.add("SRR6233828"); testCases.add("SRR6233829");
        //String command = "macs2 callpeak -t /Users/yuxiangzhang/Desktop/\"Research paper\"/";
        String command = "macs2 callpeak -t /Users/yuxiangzhang/Desktop/Pc3/";
        for (String commands: testCases){
            System.out.println(command + commands + "/" + commands + ".sam --nomodel --extsize 800 -f SAM -g 2.9e9 --outdir " + commands + " -n " + commands);
        }
    }
    static void peakCallPairedEnd(){
        ArrayList<String>testCases = new ArrayList<String>();
        testCases.add("SRR5018853"); testCases.add("SRR7968966"); testCases.add("SRR7968967");
        //String command = "macs2 callpeak -t /Users/yuxiangzhang/Desktop/\"Research paper\"/";
        String command = "macs2 callpeak -t /Users/yuxiangzhang/Desktop/\"Research paper\"/\"Paired End\"/";
        for (String commands: testCases){
            System.out.println(command + commands + "/" + commands + ".sam --nomodel --extsize 800 -f SAM -g 2.9e9 --outdir " + commands + " -n " + commands);
        }
    }
    static void convertToCorrectFormat()throws Exception{
        Scanner fileRead = new Scanner(new FileReader("sequences.txt"));
        ArrayList<ArrayList<String>> list = new ArrayList<>();
        while (fileRead.hasNext()) {
            ArrayList<String>line = new ArrayList<>();
            //StringTokenizer st = new StringTokenizer(fileRead.nextLine());

            String[]arr = fileRead.nextLine().split("	");
            line.add(arr[0]);
            line.add(arr[2].substring(0, arr[2].indexOf(",")));

            //System.out.println(Arrays.toString(arr) + " Size: " + arr.length);
            list.add(line);
        }

        ArrayList<String>testCases = new ArrayList<>();
        /*
        testCases.add("SRR6493909"); testCases.add("SRR6493910"); testCases.add("SRR5280434"); testCases.add("SRR5280435");
        testCases.add("SRR7968891"); testCases.add("SRR7968892"); testCases.add("SRR7192349");testCases.add("SRR10820724");
        testCases.add("SRR10820725"); testCases.add("SRR6233822"); testCases.add("SRR6233823"); testCases.add("SRR6233824"); testCases.add("SRR6233825");
        testCases.add("SRR6233826"); testCases.add("SRR6233827"); testCases.add("SRR6233828"); testCases.add("SRR6233829");
         */
        testCases.add("Huvec"); testCases.add("SRR7968966"); testCases.add("SRR7968967");

        String prefix = "/Users/yuxiangzhang/Desktop/Research paper/";
        for(String command: testCases){
            String fileName = prefix + command + "/" + command + "/" + command + "_peaks.narrowPeak";
            System.out.println(fileName);
            Scanner in = new Scanner(new FileReader(fileName));

            PrintWriter out = new PrintWriter(command + "_correctFormat.narrowPeak");

            while(in.hasNext()) {
                boolean possible = false;
                String s = in.nextLine();
                StringTokenizer st = new StringTokenizer(s);
                String query = st.nextToken();
                for (int i = 0; i < list.size(); i ++) {
                    ArrayList<String>newLine = list.get(i);
                    //System.out.println(query + " " + newLine + " " + newLine.contains(query));
                    if(newLine.contains(query)) {
                        out.println(s.replaceAll(query, newLine.get(0)));
                        possible = true;
                        break;
                    }
                }

                if(!possible) System.out.println(query);
                out.flush();
            }

            out.close();
        }
    }
    static void check() throws FileNotFoundException {
        String path1 = "/Users/yuxiangzhang/Desktop/Research paper/501A/test.bed";
        String path2 = "/Users/yuxiangzhang/Desktop/Research paper/501A/501Afinal.bed";
        Scanner file1 = new Scanner(new FileReader(path1));
        Scanner file2 = new Scanner(new FileReader(path2));
        HashSet<String>set = new HashSet<>();
        while(file1.hasNext()){
            String s = file1.nextLine().replaceAll("\t", "");
            set.add(s);
            //System.out.println(s);
        }

        while(file2.hasNext()){
            String s = file2.nextLine().replaceAll(" ", "");
            if(!set.contains(s)) {
                System.out.println(s);
            }
            else System.out.println("correct");
        }

        file1.close();
        file2.close();
    }
    static void prune() throws Exception {
        Scanner pathIn = new Scanner(new FileReader("sequences.txt"));
        HashMap<String, Long> map = new HashMap<>();
        while (pathIn.hasNext()) {
            String[] arr = pathIn.nextLine().split("\t");
            //System.out.println(Arrays.toString(arr));
            map.put(arr[0], Long.parseLong(arr[1].replaceAll(",", "").trim()));
        }

        ArrayList<String> paths = new ArrayList<String>();
        /*
        paths.add("501A");
        paths.add("Breast");
        paths.add("HKC8");
        paths.add("K562");
        paths.add("SRR7192349");
        paths.add("Pc3");
        */
         paths.add("HepG2"); paths.add("Huvec");
        String prefix = "/Users/yuxiangzhang/Desktop/Research paper/";
        String suffix = "sormer.bed";
        for (String s1 : paths) {
            String path = prefix + s1 + "/" + s1 + suffix;
            Scanner in = new Scanner(new FileReader(path));

            long count = 1;
            PrintWriter out = new PrintWriter(prefix + s1 + "final.bed");
            while (in.hasNext()) {
                StringTokenizer st = new StringTokenizer(in.nextLine());
                String s = st.nextToken();
                Long a = Long.parseLong(st.nextToken());
                Long b = Long.parseLong(st.nextToken());
                if (a > map.get(s) || b > map.get(s)) {
                    System.out.println(s + " " + a + " " + b);
                    continue;
                }
                if ((b - a) <= 75000) {
                    out.println(s + "   " + a + "   " + b + "   " + s1 + "_peak_" + count);
                    count++;
                }
                else System.out.println(a + " " + b);
                out.flush();
            }

            out.close();
        }
    }

}