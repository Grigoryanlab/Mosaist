#include "msttypes.h"

using namespace std;
using namespace MST;

void tokenEqualityFailure(string testName, int tokenInd, string expectedToken, string observedToken) {
  MstUtils::error("failed on '" + testName + "': token " + to_string(tokenInd) + " was expected to be <STRING>" + expectedToken + "</STRING> but was actually <STRING>" + observedToken + "</STRING>");
}

void tokenCardinalityFailure(string testName, int numExpectedTokens, int numObservedTokens) {
  MstUtils::error("failed on '" + testName + "': expected to find " + to_string(numExpectedTokens) + " tokens but actually found " + to_string(numObservedTokens));
}

void thrownExceptionFailure(string testName) {
  MstUtils::error("failed on '" + testName + "': expected to return a string but actually threw an exception");
}

int main(int argc, char** argv) {
  map<string, pair<string, vector<string>>> testStrings;
  testStrings["no quotes"] = make_pair(
    "I am a simple string with no quotes",
    vector<string>({"I", "am", "a", "simple", "string", "with", "no", "quotes"})
  );
  testStrings["quote at start"] = make_pair(
    "'I am a simple string with no quotes' - a string",
    vector<string>({"I am a simple string with no quotes", "-", "a", "string"})
  );
  testStrings["quote at end"] = make_pair(
    "A string once said, 'I am a simple string with no quotes'",
    vector<string>({"A", "string", "once", "said,", "I am a simple string with no quotes"})
  );
  testStrings["quote in middle"] = make_pair(
    "To quote my stringy friend, 'I am a simple string with no quotes', although that is not the case for myself",
    vector<string>({"To", "quote", "my", "stringy", "friend,", "I am a simple string with no quotes,", "although", "that", "is", "not", "the", "case", "for", "myself"})
  );
  /*testStrings["contraction"] = make_pair(
    "Don't include contractions if you want quote awareness",
    vector<string>()
  );*/
  testStrings["double quote"] = make_pair(
    "\"Double quotes are quotes too!\"",
    vector<string>({"Double quotes are quotes too!"})
  );
  testStrings["command quote"] = make_pair(
    "Commands are usually quoted like `command arg1 arg2 ...` to differentiate them from other quotations",
    vector<string>({"Commands", "are", "usually", "quoted", "like", "command arg1 arg2 ...", "to", "differentiate", "them", "from", "other", "quotations"})
  );
  testStrings["mixed quotes"] = make_pair(
    "Commands can include quotations of their own: `command 'arg1 with spaces' arg2 ...` and this should work fine",
    vector<string>({"Commands", "can", "include", "quotations", "of", "their", "own:", "command 'arg1 with spaces' arg2 ...", "and", "this", "should", "work", "fine"})
  );
  testStrings["contraction in quote"] = make_pair(
    "Some quotations, like \"don't\", have contractions",
    vector<string>({"Some", "quotations,", "like", "don't,", "have", "contractions"})
  );
  testStrings["only quote"] = make_pair(
    "'This is an unattributed quotation'",
    vector<string>({"This is an unattributed quotation"})
  );
  testStrings["empty quote"] = make_pair(
    "\"\"",
    vector<string>({""})
  );
  testStrings["escaped quote"] = make_pair(
    "Sometimes you want a literal quote like \\\" or many like \\\' \\\" \\\' and this should work fine",
    vector<string>({"Sometimes", "you", "want", "a", "literal", "quote", "like", "\"", "or", "many", "like", "'", "\"", "'", "and", "this", "should", "work", "fine"})
  );
  testStrings["quote with trailing characters"] = make_pair(
    "Like in Bash, characters before and after quotes with no intervening delimiter are concatenated, making a'bc'd equivalent to abcd",
    vector<string>({"Like", "in", "Bash,", "characters", "before", "and", "after", "quotes", "with", "no", "intervening", "delimiter", "are", "concatenated,", "making", "abcd", "equivalent", "to", "abcd"})
  );
  testStrings["quote with escaped backslashes"] = make_pair(
    "Backslashes can be escaped like so: \\\\",
    vector<string>({"Backslashes", "can", "be", "escaped", "like", "so:", "\\"})
  );
  testStrings["escaping with MstUtils::escape"] = make_pair(
    MstUtils::escape("Let's make sure MstUtils::escape works correctly", "'"),
    vector<string>({"Let's", "make", "sure", "MstUtils::escape", "works", "correctly"})
  );
  testStrings["no delimiters"] = make_pair(
    "No delim? Is 'abc' grouped and \\\\ unescaped?",
    vector<string>({"N", "o", " ", "d", "e", "l", "i", "m", "?", " ", "I", "s", " ", "abc", " ", "g", "r", "o", "u", "p", "e", "d", " ", "a", "n", "d", " ", "\\", " ", "u", "n", "e", "s", "c", "a", "p", "e", "d", "?"})
  );

  string quoteMarks = "'\"`";
  string delimiters = " ";
  bool skipTrailingDelims = true;
  for (auto i = testStrings.begin(); i != testStrings.end(); ++i) {
    string testName = i->first, testString = i->second.first;
    vector<string> testTokens = i->second.second;
    string testDelims = MstUtils::stringsEqual(testName, "no delimiters") ? "" : delimiters;
    cout << "Testing '" + testName + "': " + testString << endl;
    int j = 0;
    try {
      while (!testString.empty()) {
        char firstChar = testString[0];
        string token = MstUtils::nextQuoteAwareToken(testString, quoteMarks, testDelims, skipTrailingDelims);
        if (!MstUtils::stringsEqual(token, testTokens[j])) tokenEqualityFailure(testName, j, testTokens[j], token);
        j++;
      }
      if (j != testTokens.size()) tokenCardinalityFailure(testName, testTokens.size(), j);
      if (MstUtils::stringsEqual(testName, "contraction")) {
        MstUtils::error("failed on 'contraction': expected to have thrown an exception but actually returned a string");
      }
    } catch (...) {
      if (!MstUtils::stringsEqual(testName, "contraction")) thrownExceptionFailure(testName);
    }
  }

  cout << "Done: all tests passed!" << endl;
  return 0;
}
