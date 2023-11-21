import 'dart:async';
import 'dart:convert';
import 'package:flutter/material.dart';
import 'package:http/http.dart' as http;
import './rotte.dart';
import './results.dart';

Future<List<Map<String, dynamic>>> richiestaDati(
    String id, String classs) async {
  final response = await http.get(Uri.parse(
      '$apiUrl_grigliadipartenza/grigliadipartenza?id=$id&class=$classs'));

  if (response.statusCode == 200) {
    // If the server did return a 200 OK response,
    // then parse the JSON.
    return List<Map<String, dynamic>>.from(jsonDecode(response.body));
  } else {
    // If the server did not return a 200 OK response,
    // then throw an exception.
    throw Exception('Failed to load classes');
  }
}

class StartingGridRoute extends StatefulWidget {
  final String raceid;
  final String category;
  const StartingGridRoute(this.raceid, this.category, {Key? key})
      : super(key: key);

  @override
  // ignore: library_private_types_in_public_api
  _StartingGridRoute createState() => _StartingGridRoute();
}

class _StartingGridRoute extends State<StartingGridRoute> {
  late Future<List<Map<String, dynamic>>> futureResults;

  @override
  void initState() {
    super.initState();
    futureResults = richiestaDati(widget.raceid, widget.category);
  }

  @override
  Widget build(BuildContext context) {
    // ignore: no_leading_underscores_for_local_identifiers
    return Scaffold(
      backgroundColor: Colors.white,
      appBar: AppBar(
        leading: const BackButton(color: Colors.white),
        title: const Text('GRIGLIA DI PARTENZA',
            style: TextStyle(color: Colors.white, fontSize: 25)),
        centerTitle: true,
        backgroundColor: Colors.blueAccent,
      ),
      body: GestureDetector(
        onPanUpdate: (details) {
          Navigator.pushReplacement(
            context,
            MaterialPageRoute(
              builder: (context) =>
                  LeaderboardRoute(widget.raceid, widget.category),
            ),
          );
        },
        child: RefreshIndicator(
          color: Colors.white,
          backgroundColor: Colors.blueAccent,
          onRefresh: () async {
            Navigator.pushReplacement(
                context,
                MaterialPageRoute(
                    builder: (BuildContext context) => super.widget));
          },
          child: Center(
            child: FutureBuilder<List<Map<String, dynamic>>>(
              future: futureResults,
              builder: (context, snapshot) {
                if (snapshot.hasData) {
                  List<Map<String, dynamic>> result = snapshot.data!;
                  return ListView.builder(
                      itemCount: result.length,
                      itemBuilder: ((context, index) => RichText(
                            text: TextSpan(
                              children: <TextSpan>[
                                TextSpan(
                                  text:
                                      // ignore: prefer_interpolation_to_compose_strings
                                      "${" -----------------------------------------------------------------------------------------------------\n\n ID: " + result[index]["id"]}\n",
                                  style: const TextStyle(color: Colors.black),
                                ),
                                TextSpan(
                                  text:
                                      // ignore: prefer_interpolation_to_compose_strings
                                      "${"${" ATLETA: " + result[index]["name"]} " + result[index]["surname"]}\n",
                                  style: const TextStyle(color: Colors.black),
                                ),
                                TextSpan(
                                  text:
                                      " PARTENZA: ${DateTime.parse(result[index]["starttime"]).toLocal().toString().substring(0, DateTime.parse(result[index]["starttime"]).toLocal().toString().length - 4)}",
                                  style: const TextStyle(color: Colors.black),
                                ),
                              ],
                            ),
                          )));
                } else if (snapshot.hasError) {
                  return Text('${snapshot.error}');
                }
                // By default, show a loading spinner.
                return const CircularProgressIndicator(
                  color: Colors.blueAccent,
                  backgroundColor: Colors.red,
                );
              },
            ),
          ),
        ),
      ),
    );
  }
}
