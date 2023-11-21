import 'dart:async';
import 'dart:convert';
import 'package:flutter/material.dart';
import 'package:http/http.dart' as http;
import './rotte.dart';
import './grigliadipartenza.dart';

Future<List<Map<String, dynamic>>> richiestaRisultati(
    String id, String classs) async {
  final response =
      await http.get(Uri.parse('$apiUrl_results/results?id=$id&class=$classs'));

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

Color getColor(String date) {
  if (DateTime.now().difference(DateTime.parse(date)).inMinutes > 4) {
    return Colors.black;
  } else {
    return Colors.green;
  }
}

Color getColorStarted(String time) {
  if (time != "0") {
    return Colors.black;
  } else {
    return Colors.red;
  }
}

class LeaderboardRoute extends StatefulWidget {
  final String raceid;
  final String category;
  const LeaderboardRoute(this.raceid, this.category, {Key? key})
      : super(key: key);

  @override
  // ignore: library_private_types_in_public_api
  _LeaderboardRouteState createState() => _LeaderboardRouteState();
}

class _LeaderboardRouteState extends State<LeaderboardRoute> {
  late Future<List<Map<String, dynamic>>> futureResults;

  @override
  void initState() {
    super.initState();
    futureResults = richiestaRisultati(widget.raceid, widget.category);
  }

  @override
  Widget build(BuildContext context) {
    // ignore: no_leading_underscores_for_local_identifiers
    return Scaffold(
      backgroundColor: Colors.white,
      appBar: AppBar(
        leading: const BackButton(color: Colors.white),
        title: const Text('CLASSIFICA',
            style: TextStyle(color: Colors.white, fontSize: 30)),
        centerTitle: true,
        backgroundColor: Colors.blueAccent,
      ),
      body: GestureDetector(
        onPanUpdate: (details) {
          Navigator.pushReplacement(
            context,
            MaterialPageRoute(
              builder: (context) =>
                  StartingGridRoute(widget.raceid, widget.category),
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
                      // ignore: prefer_interpolation_to_compose_strings
                      itemBuilder: ((context, index) => RichText(
                            text: TextSpan(
                              children: <TextSpan>[
                                TextSpan(
                                  text:
                                      // ignore: prefer_interpolation_to_compose_strings
                                      "${" -----------------------------------------------------------------------------------------------------\n\n POSIZIONE: " + result[index]["posizione"]}\n",
                                  style: const TextStyle(color: Colors.black),
                                ),
                                TextSpan(
                                  // ignore: prefer_interpolation_to_compose_strings
                                  text: "${" ID: " + result[index]["id"]}\n",
                                  style: const TextStyle(color: Colors.black),
                                ),
                                TextSpan(
                                  text:
                                      // ignore: prefer_interpolation_to_compose_strings
                                      "${" ATLETA: " + result[index]["atleta"]}\n",
                                  style: const TextStyle(color: Colors.black),
                                ),
                                TextSpan(
                                  text:
                                      " INIZIO: ${DateTime.parse(result[index]["tempoInizio"]).toLocal().toString().substring(0, DateTime.parse(result[index]["tempoInizio"]).toLocal().toString().length - 4)}\n",
                                  style: const TextStyle(color: Colors.black),
                                ),
                                TextSpan(
                                  text:
                                      " ARRIVO: ${DateTime.parse(result[index]["tempoFine"]).toLocal().toString().substring(0, DateTime.parse(result[index]["tempoFine"]).toLocal().toString().length - 4)}\n",
                                  style: TextStyle(
                                      color:
                                          getColor(result[index]["tempoFine"])),
                                ),
                                TextSpan(
                                  text:
                                      // ignore: prefer_interpolation_to_compose_strings
                                      "${" TEMPO: " + result[index]["tempo"]} s \n",
                                  style: TextStyle(
                                      color: getColorStarted(
                                          result[index]["tempo"])),
                                ),
                                TextSpan(
                                    // ignore: prefer_interpolation_to_compose_strings
                                    text: " ORGANIZZAZIONE: " +
                                        result[index]["organizzazione"],
                                    style:
                                        const TextStyle(color: Colors.black)),
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
