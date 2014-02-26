package ai;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class test {

	public static void main(String[] args){
		
		String str = ")**^S99a8";
		Pattern pattern = Pattern.compile("[A-Za-z]");
		Matcher matcher = pattern.matcher(str);
		
		while(matcher.find()){
			System.out.print(matcher.start()+ " ");
			System.out.println(matcher.end());
		}
		
	}
}
