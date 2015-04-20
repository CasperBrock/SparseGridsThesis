package testing;

import grid.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.nio.file.FileSystem;
import java.util.LinkedList;
import java.util.List;

public class Experiments {

	public static void main(String[] args) {
		List<String> data = new LinkedList<String>();
		data.add(getInfo());
		writeToFile("test", data);
	}
	
	private static String getInfo() {
		RuntimeMXBean runtimeMxBean = ManagementFactory.getRuntimeMXBean();
		List<String> arguments = runtimeMxBean.getInputArguments();
		StringBuilder sb = new StringBuilder();
		sb.append("# OS: " +
		System.getProperty("os.name") + " " +
		System.getProperty("os.version") + " " +
		System.getProperty("os.arch"));
		sb.append('\n');
		sb.append("# JVM: " +
		System.getProperty("java.vendor") + " " +
		System.getProperty("java.version"));
		sb.append('\n');
		// This line works only on MS Windows:
		sb.append("# CPU: " + System.getenv("PROCESSOR_IDENTIFIER"));
		sb.append('\n');
		sb.append("# Available processors: " + Runtime.getRuntime().availableProcessors());
		sb.append('\n');
		sb.append("# Maximum memory: " + Runtime.getRuntime().maxMemory() / (1024 * 1024) + " MB");
		sb.append('\n');
		java.util.Date now = new java.util.Date();
		sb.append("# Date: " +
		new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ssZ").format(now));
		sb.append('\n');
		sb.append("# VM Arguments: ");
		for(String s : arguments)
			sb.append(s + " ");
		sb.append('\n');
		return sb.toString();
	}
	
	private static void writeToFile(String filename, List<String> data)
	{
		try {
			File dir = new File("Experiments");
			if(!dir.exists())
				dir.mkdir();
			File file = new File("Experiments" + File.separator + filename + ".dat");
 
			// if file exists deletes it then creates it again
			if (file.exists()) {
				file.delete();
			}
			file.createNewFile();
			
 
			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			for(String s : data) {
				bw.write(s);
				bw.write("\n");
			}
			bw.close();
 
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
