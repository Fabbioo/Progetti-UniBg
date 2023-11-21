import 'dart:async';
import 'dart:convert';
import 'package:flutter/material.dart';
import 'package:http/http.dart' as http;
import './rotte.dart';
import './grigliadipartenza.dart';

Future<List<String>> fetchClasses(String id) async {
  final response =
      await http.get(Uri.parse('$apiUrl_listclasses/listclasses?id=$id'));

  if (response.statusCode == 200) {
    // If the server did return a 200 OK response,
    // then parse the JSON.
    return List<String>.from(jsonDecode(response.body));
  } else {
    // If the server did not return a 200 OK response,
    // then throw an exception.
    throw Exception('Failed to load classes');
  }
}

class ClassesRoute extends StatefulWidget {
  final String raceid;
  const ClassesRoute(this.raceid, {Key? key}) : super(key: key);

  @override
  // ignore: library_private_types_in_public_api
  _ClassesRouteState createState() => _ClassesRouteState();
}

class _ClassesRouteState extends State<ClassesRoute> {
  late Future<List<String>> futureClasses;

  @override
  void initState() {
    super.initState();
    futureClasses = fetchClasses(widget.raceid);
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      backgroundColor: Colors.white,
      appBar: AppBar(
        leading: const BackButton(color: Colors.white),
        title: const Text('CATEGORIE',
            style: TextStyle(color: Colors.white, fontSize: 30)),
        centerTitle: true,
        backgroundColor: Colors.blueAccent,
      ),
      body: RefreshIndicator(
        color: Colors.white,
        backgroundColor: Colors.blueAccent,
        onRefresh: () async {
          Navigator.pushReplacement(
              context,
              MaterialPageRoute(
                  builder: (BuildContext context) => super.widget));
        },
        child: Center(
          child: FutureBuilder<List<String>>(
            future: futureClasses,
            builder: (context, snapshot) {
              if (snapshot.hasData) {
                List<String> classes = snapshot.data!;
                return ListView.builder(
                    itemCount: classes.length,
                    itemBuilder: ((context, index) => ElevatedButton(
                          onPressed: () {
                            Navigator.push(
                              context,
                              MaterialPageRoute(
                                builder: (context) => StartingGridRoute(
                                    widget.raceid, classes[index]),
                              ),
                            );
                          },
                          style: ElevatedButton.styleFrom(
                            primary: Colors.blue,
                            side: const BorderSide(
                                width: 2.0, color: Colors.white),
                            shape: RoundedRectangleBorder(
                              borderRadius: BorderRadius.circular(20),
                            ),
                            padding: const EdgeInsets.symmetric(
                                horizontal: 50, vertical: 8),
                          ),
                          child: Text(classes[index],
                              style: const TextStyle(
                                  fontSize: 20, color: Colors.white)),
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
    );
  }
}
