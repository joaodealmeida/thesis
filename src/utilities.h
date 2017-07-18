#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <cctype>
using namespace std;

bool to_bool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}

void tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ") {
	// skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);

	// add empty string tokens between contiguous delimiters
	for (unsigned int i = 0; i < lastPos && string::npos != lastPos; i++)
		tokens.push_back("");

	// find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos) {
		// found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));

		// skip delimiters. Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);

		// add empty string tokens between contiguous delimiters
		for (unsigned int i = 0; i < lastPos - pos - 1 && string::npos != lastPos; i++)
			tokens.push_back("");

		// find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

vector<string> splitTSV(const string& str) {
	vector<string> tokens;

	tokenize(str, tokens, "\t");

	return tokens;
}

vector<string> splitOverlay(const string& str){
	vector<string> tokens;
	tokenize(str, tokens, " ");
	return tokens;
}

vector <string> splitCyto(const string& str){
	vector<string> tokens;
	tokenize(str, tokens, ",");
	return tokens;
}

vector <string> splitPathway(const string& str){
	vector<string> tokens;
	tokenize(str, tokens, "|");
	return tokens;
}

vector <string> splitDoubleComma(const string& str){
	vector<string> tokens;
	tokenize(str, tokens, ";");
	return tokens;
}


bool isPresent(const vector<string>& v, const string& geneName){
	return (std::find(v.begin(), v.end(), geneName) != v.end());
}
/*
* String manipulation functions
*/

string str_tolower(string str) {
	transform(str.begin(), str.end(), str.begin(), ::tolower);

	return str;
}

string str_trimChar(string str, char ch) {
	for (unsigned int i = 0; i < str.size(); i++) {
		while (str[i] == ch)
			str.erase(str.begin() + i);
	}

	return str;
}

string str_replaceChar(string str, char ch, char newCh) {
	for (unsigned int i = 0; i < str.size(); i++) {
		if (str[i] == ch)
			str[i] = newCh;
	}

	return str;
}

string str_trimEqualAdjacentChars(string str, char ch) {
	for (unsigned int i = 1; i < str.size(); i++) {
		while (str[i - 1] == ch && str[i] == ch)
			str.erase(str.begin() + i);
	}

	return str;
}

std::string ReplaceAll(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

string str_slug(string str) {
	str = str_tolower(str);

	str = str_trimChar(str, '(');
	str = str_trimChar(str, ')');

	str = str_replaceChar(str, ' ', '-');

	str = str_trimEqualAdjacentChars(str, '-');

	return str;
}
